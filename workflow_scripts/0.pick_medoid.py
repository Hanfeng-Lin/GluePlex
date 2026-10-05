#!/usr/bin/env python
"""Pick which experimental chain to dock, for the `medoid` and `pesto` arms.

Those two arms union PeSTo actives over every structure of the target and then
dock ONE of them. Because the actives are already unioned, the docked chain no
longer affects the restraints at all -- it only sets the conformer. This script
makes that choice reproducible instead of leaving it to the eye.

  medoid   the chain with the lowest mean CA RMSD to all the others: the
           conformer that represents the ensemble rather than any extreme of
           it, and no scoring model is involved.

  spread   the largest pairwise RMSD in the set. Read it before trusting the
           pick: when every candidate is within ~0.5 A of every other the
           choice is nearly inconsequential, and only a set with real
           conformational range makes it worth arguing about.

Chains are compared over the residues they have in common. Residues are paired
by NUMBER whenever the numbering agrees, which for chains of the same protein
cut to the same range is exact; sequence alignment is used only when it does
not, because it is strictly worse when numbering does agree -- a global
alignment can slide a gap a few positions and match two residues ten apart in
space, and those few atoms then dominate the RMSD.

A chain covering far fewer residues than the rest is excluded by default
(`--min-coverage`): it would otherwise look central simply by having little to
disagree about.

    python workflow_scripts/0.pick_medoid.py structures/*.pdb
    python workflow_scripts/0.pick_medoid.py structures/*.pdb --emit-config
"""
import argparse
import itertools
import os
import sys

import numpy as np
from Bio.PDB import PDBParser, Superimposer
from Bio.Align import PairwiseAligner
from Bio.Data.IUPACData import protein_letters_3to1

_QUIET = PDBParser(QUIET=True)


def ca_atoms(path, chain_id=None):
    """(CA atoms, one-letter sequence, chain id) for one chain of a PDB."""
    st = _QUIET.get_structure("s", path)
    model = st[0]
    if chain_id is None:
        chains = [c for c in model if any(r.get_id()[0] == " " for r in c)]
        if not chains:
            return [], "", None
        # the longest chain, so a PDB carrying a small partner still gives the
        # target rather than whichever chain happens to be first
        chain = max(chains, key=lambda c: sum(1 for r in c if r.get_id()[0] == " "))
    else:
        if chain_id not in model:
            return [], "", chain_id
        chain = model[chain_id]
    atoms, seq = [], []
    for res in chain.get_residues():
        if res.get_id()[0] != " ":
            continue
        ca = res.child_dict.get("CA")
        if ca:
            atoms.append(ca)
            seq.append(protein_letters_3to1.get(res.get_resname().title(), "X"))
    return atoms, "".join(seq), chain.get_id()


def pair_by_sequence(a1, s1, a2, s2):
    """Match two chains residue-by-residue through a global sequence alignment."""
    al = PairwiseAligner()
    al.mode = "global"
    al.match_score, al.mismatch_score = 2, -1
    al.open_gap_score, al.extend_gap_score = -5, -0.5
    aln = al.align(s1, s2)[0]
    p1, p2, i, j = [], [], 0, 0
    for c1, c2 in zip(aln[0], aln[1]):
        if c1 != "-" and c2 != "-":
            p1.append(a1[i])
            p2.append(a2[j])
        if c1 != "-":
            i += 1
        if c2 != "-":
            j += 1
    return p1, p2


def pair_by_number(a1, a2):
    """Match two chains on author residue numbering."""
    n1 = {at.get_parent().get_id()[1]: at for at in a1}
    n2 = {at.get_parent().get_id()[1]: at for at in a2}
    common = sorted(set(n1) & set(n2))
    return [n1[r] for r in common], [n2[r] for r in common]


def pair_atoms(a1, s1, a2, s2, number_frac=0.8):
    """Pair by residue number if the numbering lines up, else by sequence.

    Measured on the 12 ALK chains this was written for: numbering agrees on
    214-216 of ~220 residues, and pairing on it reproduces the pairwise RMSDs
    exactly. Sequence alignment on the same chains misplaces two residues per
    pair by ten positions and inflates one RMSD from 0.49 to 1.85 A -- enough
    to move which chain looks most central.
    """
    p1, p2 = pair_by_number(a1, a2)
    if len(p1) >= number_frac * min(len(a1), len(a2)):
        return p1, p2
    return pair_by_sequence(a1, s1, a2, s2)


def rmsd(a1, s1, a2, s2, min_pairs=10):
    p1, p2 = pair_atoms(a1, s1, a2, s2)
    if len(p1) < min_pairs:
        return float("nan")
    sup = Superimposer()
    sup.set_atoms(p1, p2)
    return float(sup.rms)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pdbs", nargs="+", help="candidate structures of ONE target")
    ap.add_argument("--chain", default=None,
                    help="chain id to use in every file; default = the longest")
    ap.add_argument("--min-coverage", type=float, default=0.9,
                    help="drop chains with fewer than this fraction of the "
                         "largest residue count (default 0.9)")
    ap.add_argument("--emit-config", action="store_true",
                    help="print the target_pdb / pesto_extra_pdbs block")
    a = ap.parse_args()

    loaded = []
    for p in a.pdbs:
        atoms, seq, cid = ca_atoms(p, a.chain)
        if len(atoms) < 10:
            print(f"  skipping {os.path.basename(p)}: "
                  f"{len(atoms)} CA atoms", file=sys.stderr)
            continue
        loaded.append({"path": p, "atoms": atoms, "seq": seq,
                       "chain": cid, "n": len(atoms)})
    if len(loaded) < 2:
        sys.exit("Need at least 2 usable structures to pick a medoid.")

    nmax = max(x["n"] for x in loaded)
    kept = [x for x in loaded if x["n"] >= a.min_coverage * nmax]
    for x in loaded:
        if x not in kept:
            print(f"  excluded {os.path.basename(x['path'])}: {x['n']}/{nmax} "
                  f"residues, below --min-coverage {a.min_coverage}",
                  file=sys.stderr)
    if len(kept) < 2:
        sys.exit("Fewer than 2 structures survived --min-coverage.")

    n = len(kept)
    D = np.zeros((n, n))
    for i, j in itertools.combinations(range(n), 2):
        r = rmsd(kept[i]["atoms"], kept[i]["seq"],
                 kept[j]["atoms"], kept[j]["seq"])
        D[i, j] = D[j, i] = r
    if np.isnan(D).any():
        print("  warning: some pairs had too few aligned residues to compare",
              file=sys.stderr)

    mean = np.nanmean(np.where(np.eye(n, dtype=bool), np.nan, D), axis=1)
    best = int(np.nanargmin(mean))
    spread = float(np.nanmax(D))

    print(f"\n{n} structures, pairwise CA RMSD over matched residues\n")
    print(f"  {'structure':34}{'chain':>6}{'res':>6}{'mean RMSD':>11}")
    print("  " + "-" * 57)
    for i in np.argsort(mean):
        mark = "  <- medoid" if i == best else ""
        print(f"  {os.path.basename(kept[i]['path']):34}"
              f"{str(kept[i]['chain']):>6}{kept[i]['n']:>6}"
              f"{mean[i]:>11.3f}{mark}")
    print(f"\n  spread (largest pairwise RMSD): {spread:.3f} A")
    if spread < 0.5:
        print("  -> the candidates are nearly identical; which one you dock "
              "should not matter much.")
    else:
        print("  -> a real conformational range, so the choice is worth making "
              "deliberately.")

    if a.emit_config:
        med = kept[best]["path"]
        print("\n# paste into your arm_medoid.yaml")
        print(f'target_pdb: "{med}"')
        print("pesto_extra_pdbs:")
        for x in kept:
            if x["path"] != med:
                print(f'  - "{x["path"]}"')
        print("pesto_aggregate: union")


if __name__ == "__main__":
    main()
