"""Prepare pass-2 Boltz jobs from top-N step-3 predictions.

Reads the step-3 atom.csv, sorts by a rank metric, selects the top-N rows,
converts each selected PDB to mmCIF via ChimeraX (Boltz's PDB template parser
fails on bare-ATOM PDBs without SEQRES), and emits one Boltz input YAML +
one wrapper Python run script per selection into the step4 working directory.
Does NOT execute Boltz — see 4.2.run_second_pass.py.
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import yaml


# Wrapper script template: mirrors 3.1.run_batch_clusters.SCRIPT_TEMPLATE so
# the pass-2 jobs run through the same `boltz predict` CLI the rest of the
# pipeline uses. Python {{ }} is escaped for .format().
SCRIPT_TEMPLATE = """import os
import subprocess
import sys

YAML_CONTENT = \"\"\"
version: 1
sequences:
  - protein:
      id: A  # {protein_a_name}
      sequence: {protein_a_seq}
      msa: empty
  - protein:
      id: B  # {protein_b_name}
      sequence: {protein_b_seq}
      msa: empty
  - ligand:
      id: C
      smiles: '{ligand_smiles}'
templates:
    - cif: {template_cif}
      chain_id: [A, B]
      template_id: [A, B]
\"\"\"

def main():
    # the launcher may pin the device through the environment; only fall back
    # to the value baked in at generation time if it did not
    os.environ.setdefault("CUDA_VISIBLE_DEVICES", str({gpu_id}))
    yaml_file = "{yaml_filename}"
    with open(yaml_file, "w") as f:
        f.write(YAML_CONTENT)
    print(f"Created {{yaml_file}}")

    cmd = [
        "boltz",
        "predict",
        yaml_file,
        "--out_dir", "{out_dir}",
        "--output_format", "pdb",
        "--diffusion_samples", str({diffusion_samples}),
    ]
    # Same optional seed as pass 1; absent, boltz draws its own as before.
    if {boltz_seed} is not None:
        cmd += ["--seed", str({boltz_seed})]
    if {override}:
        cmd.append("--override")

    print("Running:", " ".join(cmd))
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running boltz: {{e}}")
        sys.exit(1)


if __name__ == "__main__":
    main()
"""


def select_top_n(atom_csv: Path, step3_results_dir: Path, top_n: int,
                 rank_metric: str = "pair_chains_ipTM_avg") -> list:
    """Return a list of (rank, absolute_pdb_path, score) for the top-N models."""
    df = pd.read_csv(atom_csv)
    if rank_metric not in df.columns:
        raise KeyError(
            f"Rank metric '{rank_metric}' not found in {atom_csv}. "
            f"Available columns: {list(df.columns)[:10]}..."
        )
    # Keep only pass-1 rows if pass column exists (don't re-template pass-2 outputs).
    if "pass" in df.columns:
        df = df[df["pass"] == 1]
    df_sorted = df.sort_values(by=rank_metric, ascending=False).head(top_n)
    out = []
    for rank, (_, row) in enumerate(df_sorted.iterrows(), start=1):
        pdb_path = (step3_results_dir / row["file"]).resolve()
        if not pdb_path.exists():
            print(f"Warning: top-{rank} PDB not found at {pdb_path}", file=sys.stderr)
            continue
        out.append((rank, pdb_path, float(row[rank_metric])))
    return out


def convert_pdbs_to_cif(selection, out_dir: Path) -> dict:
    """Convert each selected PDB to mmCIF using ChimeraX.

    Boltz's PDB template loader calls setup_entities() then converts to mmCIF
    internally, which fails when the input PDB lacks SEQRES/entity metadata
    (as with Boltz's own output). ChimeraX with `bestGuess true` produces
    mmCIF with proper polymer sequences so Boltz's parser succeeds.

    Returns dict: {rank -> absolute CIF path}.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    cxc_lines = []
    cif_paths = {}
    for rank, pdb_path, _ in selection:
        cif_path = (out_dir / f"{pdb_path.stem}.cif").resolve()
        cif_paths[rank] = cif_path
        cxc_lines.append(f"open {pdb_path}")
        cxc_lines.append(f"save {cif_path} format mmcif bestGuess true")
        cxc_lines.append("close session")
    cxc_lines.append("exit")

    cxc_file = out_dir / "_pdb_to_cif.cxc"
    cxc_file.write_text("\n".join(cxc_lines) + "\n")

    print(f"Running ChimeraX to convert {len(selection)} PDB templates to mmCIF...")
    try:
        subprocess.run(["chimerax", "--nogui", str(cxc_file)], check=True)
    except FileNotFoundError:
        print("Error: 'chimerax' executable not found in PATH. "
              "Install ChimeraX or add it to PATH — see README.", file=sys.stderr)
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"Error running ChimeraX: {e}", file=sys.stderr)
        sys.exit(1)

    # Verify all CIFs were produced
    missing = [r for r, p in cif_paths.items() if not p.exists()]
    if missing:
        print(f"Error: ChimeraX did not produce CIFs for ranks {missing}", file=sys.stderr)
        sys.exit(1)
    return cif_paths


def write_job_files(selection, cif_paths: dict, step3_config, step4_dir: Path,
                    diffusion_samples: int, gpus: list, override: bool,
                    boltz_seed=None) -> list:
    """Write per-template YAML + wrapper run script. Return the run script paths."""
    protein_a_name = step3_config["protein_a_name"]
    protein_b_name = step3_config["protein_b_name"]
    protein_a_seq = step3_config["protein_a_sequence"]
    protein_b_seq = step3_config["protein_b_sequence"]
    ligand_smiles = step3_config["ligand_smiles"]

    out_subdir = step4_dir / f"{protein_a_name}_{protein_b_name}_pass2_noMSA_{diffusion_samples}runs"

    scripts = []
    for idx, (rank, template_pdb, score) in enumerate(selection):
        gpu_id = gpus[idx % len(gpus)]
        job_tag = f"pass2_top{rank}"
        yaml_filename = f"{protein_a_name}_{protein_b_name}_{job_tag}.yaml"
        script_name = step4_dir / f"run_boltz2_{job_tag}.py"
        template_cif = cif_paths[rank]

        content = SCRIPT_TEMPLATE.format(
            boltz_seed=boltz_seed,
            protein_a_name=protein_a_name,
            protein_a_seq=protein_a_seq,
            protein_b_name=protein_b_name,
            protein_b_seq=protein_b_seq,
            ligand_smiles=ligand_smiles,
            template_cif=str(template_cif),
            yaml_filename=yaml_filename,
            out_dir=str(out_subdir),
            gpu_id=gpu_id,
            diffusion_samples=diffusion_samples,
            override=override,
        )
        script_name.write_text(content)
        scripts.append(script_name)
        print(f"[{job_tag}] rank={rank} score={score:.4f} gpu={gpu_id} "
              f"template={template_pdb.name} -> {script_name.name}")

    return scripts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="project_config.yaml")
    parser.add_argument("--step3_dir", required=True,
                        help="step3_boltz directory (contains atom.csv and the flattened results dir)")
    parser.add_argument("--step4_dir", required=True, help="step4_boltz directory to write outputs")
    parser.add_argument("--top_n", type=int, default=None,
                        help="Override step4.top_n from config.")
    parser.add_argument("--diffusion_samples", type=int, default=None,
                        help="Override step4.diffusion_samples from config.")
    parser.add_argument("--gpus", type=int, nargs="+", default=None,
                        help="Override step4.gpus from config.")
    parser.add_argument("--override", action="store_true",
                        help="Pass --override to the Boltz CLI in each run script.")
    parser.add_argument("--dry_run", action="store_true",
                        help="Print selected templates but do not write any files.")
    args = parser.parse_args()

    with open(args.config) as f:
        full_config = yaml.safe_load(f)
    step3_config = full_config.get("step3") or {}
    step4_config = full_config.get("step4") or {}
    if not step3_config:
        print("Error: 'step3' section missing in config.", file=sys.stderr)
        sys.exit(1)

    top_n = args.top_n if args.top_n is not None else step4_config.get("top_n", 10)
    diffusion_samples = (args.diffusion_samples if args.diffusion_samples is not None
                         else step4_config.get("diffusion_samples",
                                               step3_config.get("diffusion_samples", 20)))
    gpus = args.gpus or step4_config.get("gpus") or step3_config.get("gpus") or [0]
    rank_metric = step4_config.get("rank_metric", "pair_chains_ipTM_avg")

    step3_dir = Path(args.step3_dir).resolve()
    step4_dir = Path(args.step4_dir).resolve()
    atom_csv = step3_dir / "atom.csv"
    if not atom_csv.exists():
        print(f"Error: step-3 atom.csv not found at {atom_csv}", file=sys.stderr)
        sys.exit(1)

    protein_a_name = step3_config["protein_a_name"]
    protein_b_name = step3_config["protein_b_name"]
    step3_results_dir = step3_dir / f"{protein_a_name}-{protein_b_name}"
    if not step3_results_dir.exists():
        print(f"Error: step-3 results dir not found at {step3_results_dir}", file=sys.stderr)
        sys.exit(1)

    selection = select_top_n(atom_csv, step3_results_dir, top_n, rank_metric=rank_metric)
    if not selection:
        print("No top-N models selected; aborting.", file=sys.stderr)
        sys.exit(1)

    print(f"Selected top-{len(selection)} from {atom_csv} by {rank_metric}:")
    for rank, pdb, score in selection:
        print(f"  rank {rank:>2}  {rank_metric}={score:.4f}  {pdb.name}")

    if args.dry_run:
        print("Dry run — no files written.")
        return

    step4_dir.mkdir(parents=True, exist_ok=True)
    cif_paths = convert_pdbs_to_cif(selection, step4_dir / "templates_cif")
    scripts = write_job_files(
        selection, cif_paths, step3_config, step4_dir,
        diffusion_samples=int(diffusion_samples),
        gpus=list(gpus),
        override=bool(args.override),
        boltz_seed=full_config.get("boltz_seed"),
    )
    # Emit a manifest so 4.2 can find the run scripts deterministically.
    manifest = step4_dir / "pass2_manifest.txt"
    manifest.write_text("\n".join(str(s) for s in scripts) + "\n")
    print(f"Wrote {len(scripts)} run scripts; manifest at {manifest}")

    # Save rename mapping log
    protein_a_name = step3_config["protein_a_name"]
    protein_b_name = step3_config["protein_b_name"]
    log_path = step4_dir / "pass2_template_mapping.log"
    with open(log_path, "w") as lf:
        lf.write("Pass-2 template mapping (step 3 source -> step 4 rename)\n")
        lf.write(f"{'Rank':<6}{'pair_ipTM':>12}  {'Source (step 3)':<55} {'Pass-2 name'}\n")
        lf.write(f"{'-'*110}\n")
        for rank, template_pdb, score in selection:
            job_tag = f"pass2_top{rank}"
            new_name = f"{protein_a_name}_{protein_b_name}_{job_tag}"
            lf.write(f"{rank:<6}{score:>12.4f}  {template_pdb.name:<55} {new_name}\n")
    print(f"Template mapping log saved to {log_path}")


if __name__ == "__main__":
    main()
