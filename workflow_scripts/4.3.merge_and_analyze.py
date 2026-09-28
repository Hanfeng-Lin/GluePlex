"""Merge pass-1 + pass-2 predictions and run unified post-analysis.

Symlinks pass-1 PDBs/JSONs from `step3_boltz/{protein_a}-{protein_b}/` and
pass-2 predictions from `step4_boltz/.../predictions/**/` into
`step4_boltz/merged_pdbs/`, then invokes 3.2.post_analysis.py on that merged
directory to produce `atom_2ndpass.csv` + combined UMAP/PCA plots.
"""

import argparse
import glob
import os
import subprocess
import sys
from pathlib import Path

import yaml


def symlink_if_missing(src: Path, dst_dir: Path) -> None:
    dst = dst_dir / src.name
    if dst.exists() or dst.is_symlink():
        return
    os.symlink(src.resolve(), dst)


def populate_merged(step3_results_dir: Path, step4_dir: Path, merged_dir: Path) -> tuple:
    merged_dir.mkdir(parents=True, exist_ok=True)

    # Pass-1: flattened step3 CRBN-VAV1 dir already contains *.pdb and *.json.
    p1_pdbs = sorted(step3_results_dir.glob("*.pdb"))
    p1_jsons = sorted(step3_results_dir.glob("*.json"))
    for f in p1_pdbs + p1_jsons:
        symlink_if_missing(f, merged_dir)

    # Pass-2: scan step4_dir recursively for predictions/**/*.pdb + confidence_*.json
    p2_pdbs = [Path(p) for p in glob.glob(str(step4_dir / "**" / "*.pdb"), recursive=True)
               if "predictions" in Path(p).parts]
    p2_jsons = [Path(p) for p in glob.glob(str(step4_dir / "**" / "confidence_*.json"), recursive=True)
                if "predictions" in Path(p).parts]
    # Exclude merged_pdbs itself (would be recursive) — already filtered by "predictions" check.
    for f in p2_pdbs + p2_jsons:
        symlink_if_missing(f, merged_dir)

    return (len(p1_pdbs), len(p1_jsons), len(p2_pdbs), len(p2_jsons))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--step3_dir", required=True)
    parser.add_argument("--step4_dir", required=True)
    parser.add_argument("--output_csv", default="atom_2ndpass.csv")
    parser.add_argument("--ref_pdb", default=None, help="Reference PDB for RMSD calculation.")
    parser.add_argument("--dim_reduction", default="all", choices=["umap", "pca", "all"])
    args = parser.parse_args()

    with open(args.config) as f:
        full_config = yaml.safe_load(f)
    step3_config = full_config.get("step3") or {}
    protein_a = step3_config.get("protein_a_name", "A")
    protein_b = step3_config.get("protein_b_name", "B")
    ligand_name = step3_config.get("ligand_name", "Ligand")
    protein_a_range = step3_config.get("protein_a_residue_range", "")
    protein_b_range = step3_config.get("protein_b_residue_range", "")

    step3_dir = Path(args.step3_dir).resolve()
    step4_dir = Path(args.step4_dir).resolve()
    step3_results_dir = step3_dir / f"{protein_a}-{protein_b}"
    if not step3_results_dir.exists():
        print(f"Error: {step3_results_dir} not found — step 3 must complete first.",
              file=sys.stderr)
        sys.exit(1)

    merged_dir = step4_dir / "merged_pdbs"
    n1p, n1j, n2p, n2j = populate_merged(step3_results_dir, step4_dir, merged_dir)
    print(f"Merged: pass-1 {n1p} PDBs / {n1j} JSONs, pass-2 {n2p} PDBs / {n2j} JSONs "
          f"into {merged_dir}")
    if n2p == 0:
        print("Warning: no pass-2 PDBs found. Did 4.2 run?", file=sys.stderr)

    script_32 = Path(__file__).resolve().parent / "3.2.post_analysis.py"
    results_dir = step4_dir / f"{protein_a}-{protein_b}_2ndpass"
    cmd = [
        sys.executable, str(script_32),
        "--source_dir", str(merged_dir),
        "--results_dir", str(results_dir),
        "--output_csv", args.output_csv,
        "--align_chain", "A",
        "--extract_chain", "B",
        "--ligand_name", ligand_name,
        "--protein_name", protein_b,
        "--protein_a_name", protein_a,
        "--protein_a_residue_range", protein_a_range,
        "--protein_b_residue_range", protein_b_range,
        "--dim_reduction", args.dim_reduction,
    ]
    if args.ref_pdb:
        cmd.extend(["--ref_pdb", args.ref_pdb])
        # Multi-domain receptors (e.g. CRBN's Lon/CULT hinge) need the superposition
        # restricted to the binding domain, or domain motion inflates the interface RMSD.
        ref_align_range = full_config.get("ref_align_range")
        if ref_align_range:
            cmd.extend(["--align_range", str(ref_align_range)])
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True, cwd=str(step4_dir))
    print(f"Combined post-analysis complete. See {step4_dir / args.output_csv}")


if __name__ == "__main__":
    main()
