"""Compute per-model chain B CA RMSD against a reference PDB.

Aligns each model's chain A to the reference chain A (sequence alignment +
BioPython Superimposer), then measures chain B CA displacement without
re-superposing chain B.  Adds a `ref_rmsd` column to the input CSV in place.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, Superimposer
from Bio.Align import PairwiseAligner
from Bio.Data.IUPACData import protein_letters_3to1


def _get_ca_atoms(structure, chain_id):
    atoms, seq = [], []
    for res in structure[0][chain_id].get_residues():
        if res.get_id()[0] == " ":
            ca = res.child_dict.get("CA")
            if ca:
                atoms.append(ca)
                seq.append(protein_letters_3to1.get(res.get_resname().title(), "X"))
    return atoms, "".join(seq)


def _pair_by_sequence(atoms1, seq1, atoms2, seq2):
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5
    aln = aligner.align(seq1, seq2)[0]
    paired1, paired2 = [], []
    i, j = 0, 0
    for c1, c2 in zip(aln[0], aln[1]):
        if c1 != "-" and c2 != "-":
            paired1.append(atoms1[i])
            paired2.append(atoms2[j])
        if c1 != "-":
            i += 1
        if c2 != "-":
            j += 1
    return paired1, paired2


def compute_ref_rmsd(pdb_path, bp_parser,
                     ref_A_atoms, ref_A_seq,
                     ref_B_atoms, ref_B_seq,
                     align_chain, extract_chain):
    try:
        mod = bp_parser.get_structure("mod", pdb_path)
    except Exception:
        return float("nan")

    mod_A, mod_seq_A = _get_ca_atoms(mod, align_chain)
    p_ref_A, p_mod_A = _pair_by_sequence(ref_A_atoms, ref_A_seq, mod_A, mod_seq_A)
    if len(p_ref_A) < 10:
        return float("nan")

    sup = Superimposer()
    sup.set_atoms(p_ref_A, p_mod_A)
    sup.apply(list(mod.get_atoms()))

    mod_B, mod_seq_B = _get_ca_atoms(mod, extract_chain)
    p_ref_B, p_mod_B = _pair_by_sequence(ref_B_atoms, ref_B_seq, mod_B, mod_seq_B)
    if not p_ref_B:
        return float("nan")

    ref_coords = np.array([a.get_vector().get_array() for a in p_ref_B])
    mod_coords = np.array([a.get_vector().get_array() for a in p_mod_B])
    return float(np.sqrt(np.mean(np.sum((ref_coords - mod_coords) ** 2, axis=1))))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ref_pdb", required=True, help="Reference PDB path.")
    parser.add_argument("--input_csv", required=True, help="atom.csv to update.")
    parser.add_argument("--results_dir", required=True,
                        help="Directory containing model PDB files referenced in the CSV.")
    parser.add_argument("--align_chain", default="A")
    parser.add_argument("--extract_chain", default="B")
    args = parser.parse_args()

    bp_parser = PDBParser(QUIET=True)
    ref_st = bp_parser.get_structure("ref", args.ref_pdb)
    ref_A_atoms, ref_A_seq = _get_ca_atoms(ref_st, args.align_chain)
    ref_B_atoms, ref_B_seq = _get_ca_atoms(ref_st, args.extract_chain)
    print(f"Reference: {args.ref_pdb}")
    print(f"  Chain {args.align_chain}: {len(ref_A_atoms)} CAs")
    print(f"  Chain {args.extract_chain}: {len(ref_B_atoms)} CAs")

    df = pd.read_csv(args.input_csv)
    rmsds = []
    for _, row in df.iterrows():
        pdb_path = os.path.join(args.results_dir, row["file"])
        if not os.path.exists(pdb_path):
            rmsds.append(float("nan"))
            continue
        rmsd = compute_ref_rmsd(
            pdb_path, bp_parser,
            ref_A_atoms, ref_A_seq,
            ref_B_atoms, ref_B_seq,
            args.align_chain, args.extract_chain,
        )
        rmsds.append(rmsd)

    df["ref_rmsd"] = rmsds
    df.to_csv(args.input_csv, index=False)

    valid = [r for r in rmsds if not np.isnan(r)]
    if valid:
        print(f"ref_rmsd: {len(valid)} models, "
              f"min={min(valid):.2f}, max={max(valid):.2f}, mean={np.mean(valid):.2f}")
    else:
        print("Warning: no valid ref_rmsd computed.")
    print(f"Updated {args.input_csv} with ref_rmsd column.")


if __name__ == "__main__":
    main()
