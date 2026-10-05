#!/usr/bin/env python
"""Arm `msa`: off-the-shelf Boltz-2, no HADDOCK and no template.

Two sequences plus the SMILES, MSA from the MMseqs2 server, nothing else. This
is the condition published cofolding benchmarks measure, and it is the only
baseline that says whether the PeSTo -> HADDOCK -> template machinery buys
anything at all.

The main pipeline runs Boltz with `msa: empty` because the docking template
already supplies the fold. Dropping the template without putting the MSA back
would leave Boltz nothing to work from and would not be a fair comparison, so
the MSA goes back in here.

Reads the same config as run_glueplex.py; everything outside `step3` is ignored
because no docking happens.

    python run_msa_baseline.py --config my_target.yaml [--samples 440]
"""
import argparse
import glob
import os
import shutil
import subprocess
import sys

import yaml

# Boltz draws every diffusion sample in one batch, so a large request OOMs and
# the whole run is wasted. 110 was measured safe on a 695-residue system (40 GB
# peak on an H100). On OOM the chunk is retried at half the size, so a system
# larger than anything measured here degrades instead of failing outright.
CHUNK = 110
MIN_CHUNK = 10


def models_present(chunk_dir, tag):
    return bool(glob.glob(
        f"{chunk_dir}/boltz_results_{tag}/predictions/{tag}/*_model_*.pdb"))


def write_boltz_yaml(cfg, path):
    """Sequences + ligand, with no `templates:` key and no `msa: empty`.

    Leaving both out is what makes this the no-template arm: Boltz fills the
    MSA from the server because --use_msa_server is on the command line, and
    has no structural input at all.
    """
    s3 = cfg["step3"]
    for k in ("protein_a_sequence", "protein_b_sequence", "ligand_smiles"):
        if not s3.get(k):
            sys.exit(f"Error: step3.{k} is required for the msa arm.")
    doc = {
        "version": 1,
        "sequences": [
            {"protein": {"id": "A", "sequence": s3["protein_a_sequence"]}},
            {"protein": {"id": "B", "sequence": s3["protein_b_sequence"]}},
            {"ligand": {"id": "C", "smiles": s3["ligand_smiles"]}},
        ],
    }
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    yaml.safe_dump(doc, open(path, "w"), sort_keys=False)
    return path


def run_chunks(yml, tag, out, samples, gpu, cache, log):
    """Sample in chunks, pooling the models into one directory."""
    env = dict(os.environ)
    if gpu is not None:
        env["CUDA_VISIBLE_DEVICES"] = str(gpu)
    pool = f"{out}/models"
    os.makedirs(pool, exist_ok=True)
    boltz = shutil.which("boltz")
    if not boltz:
        sys.exit("Error: `boltz` not on PATH. Activate the boltz_analysis env.")

    n_chunks = -(-samples // CHUNK)
    for k in range(n_chunks):
        cdir = f"{out}/chunk{k:02d}"
        if models_present(cdir, tag):
            continue
        want = min(CHUNK, samples - k * CHUNK)
        while want >= MIN_CHUNK:
            shutil.rmtree(cdir, ignore_errors=True)
            cmd = [boltz, "predict", yml, "--out_dir", cdir,
                   "--output_format", "pdb",
                   "--diffusion_samples", str(want), "--use_msa_server",
                   "--msa_pairing_strategy", "greedy",
                   "--seed", str(1000 + k), "--override"]
            if cache:
                cmd += ["--cache", cache]
            mark = os.path.getsize(log) if os.path.exists(log) else 0
            with open(log, "a") as f:
                f.write(f"\n===== chunk {k} (seed {1000 + k}, n={want}) =====\n")
                f.flush()
                subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT,
                               cwd=out, env=env, text=True)
            with open(log) as f:
                f.seek(mark)
                tail = f.read()
            if "ran out of memory" in tail or not models_present(cdir, tag):
                want //= 2
                print(f"  chunk {k}: out of memory, retrying at {want} samples")
                continue
            break
        if want < MIN_CHUNK:
            sys.exit(f"Error: chunk {k} failed even at {MIN_CHUNK} samples; "
                     f"see {log}")

        src = f"{cdir}/boltz_results_{tag}/predictions/{tag}"
        for f_ in sorted(glob.glob(f"{src}/*")):
            b = os.path.basename(f_)
            if b.endswith(".npz"):
                continue
            # The chunk tag goes INSIDE the stem, never in front of it:
            # post_analysis pairs each model with confidence_<stem>.json, and a
            # leading prefix on both names breaks that lookup silently -- every
            # ipTM then reads back as 0.000 and the ranking is meaningless.
            nb = b.replace("_model_", f"_c{k:02d}_model_") if "_model_" in b \
                else f"c{k:02d}_{b}"
            shutil.move(f_, f"{pool}/{nb}")
        shutil.rmtree(cdir, ignore_errors=True)
    return pool, len(glob.glob(f"{pool}/*_model_*.pdb"))


def score(cfg, pool, out, csv, log):
    """Score with the same script the main pipeline uses.

    Keeping 3.2.post_analysis.py means the CSV carries the same columns
    (pair_chains_ipTM_avg, the chain-B CA coordinates, ref_rmsd when a
    reference is given), so the clustering and metrics read it unchanged and
    this arm stays comparable with the others.
    """
    s3 = cfg["step3"]
    # separate dirs: post_analysis copies source -> results and raises
    # SameFileError when the two are one directory
    cmd = [sys.executable, "workflow_scripts/3.2.post_analysis.py",
           "--source_dir", pool, "--results_dir", f"{out}/scored",
           "--output_csv", csv, "--align_chain", "A", "--extract_chain", "B",
           "--protein_a_name", s3.get("protein_a_name", "CRBN"),
           "--protein_name", s3.get("protein_b_name", "TARGET"),
           "--dim_reduction", "pca"]
    for key, flag in (("protein_a_residue_range", "--protein_a_residue_range"),
                      ("protein_b_residue_range", "--protein_b_residue_range")):
        if s3.get(key):
            cmd += [flag, str(s3[key])]
    # ref_pdb is optional: without it the CSV simply has no ref_rmsd column,
    # which is the normal case for a target with no solved ternary structure.
    if cfg.get("ref_pdb"):
        cmd += ["--ref_pdb", cfg["ref_pdb"]]
        if cfg.get("ref_align_range"):
            cmd += ["--align_range", str(cfg["ref_align_range"])]
    with open(log, "a") as f:
        f.write(f"\n===== scoring =====\n")
        f.flush()
        return subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT,
                              cwd=os.getcwd(), text=True).returncode


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--samples", type=int, default=440,
                    help="total diffusion samples (default matches the median "
                         "model count of a full two-pass GluePlex run)")
    ap.add_argument("--msa-cache", default=None,
                    help="Boltz --cache dir; reuses MSAs across runs of the "
                         "same sequences")
    a = ap.parse_args()

    cfg = yaml.safe_load(open(a.config))
    s3 = cfg.get("step3") or {}
    project = cfg.get("project_name", "Result_msa")
    out = os.path.join(cfg.get("base_dir", "."), project, "step_msa_notmpl")
    os.makedirs(out, exist_ok=True)
    tag = f"{s3.get('protein_b_name', 'TARGET')}_{s3.get('ligand_name', 'LIG')}"
    log = f"{out}/run.log"
    open(log, "a").close()

    gpus = s3.get("gpus") or [None]
    if len(gpus) > 1:
        print(f"Note: this runner uses one GPU; taking {gpus[0]} of {gpus}.")

    yml = write_boltz_yaml(cfg, f"{out}/{tag}.yaml")
    print(f"=== Arm msa: Boltz-2, MSA from server, no template ===")
    print(f"  {a.samples} samples in chunks of {CHUNK} -> {out}")

    pool, n = run_chunks(yml, tag, out, a.samples, gpus[0], a.msa_cache, log)
    print(f"  {n} models")
    if not n:
        sys.exit(f"Error: no models produced; see {log}")

    csv = f"{out}/atom_msa.csv"
    if score(cfg, pool, out, csv, log) or not os.path.exists(csv):
        sys.exit(f"Error: scoring failed; see {log}")
    print(f"Done. {csv}")


if __name__ == "__main__":
    main()
