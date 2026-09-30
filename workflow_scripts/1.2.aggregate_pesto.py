#!/usr/bin/env python
"""Merge PeSTo interface predictions from several structures of the same target.

Why this exists
---------------
Step 1 normally runs PeSTo on ONE structure, usually an AlphaFold model, and
everything downstream depends on the active residues it produces. That is
fragile: on the NEK7 AF2 model PeSTo's maximum score was 0.440 and exactly one
residue cleared the 0.3 threshold, leaving HADDOCK with nothing to work with.
Other experimental structures of the same domain score far higher.

This script runs the same threshold over several structures and merges the
result, writing a PeSTo-format PDB that step 2 consumes unchanged.

Two things it is careful about
------------------------------
1. **The primary structure is always one of the sources.** An aggregate built
   only from extra structures is not a superset of what the pipeline already
   had, so it can *lose* signal rather than only add to it.
2. **PeSTo renumbers its output 1..N** over observed residues, so scores from
   different structures are in mutually incompatible local numbering. Order and
   count are preserved, so each output is mapped back through its own input
   before anything is combined. Inputs must therefore share a residue-numbering
   convention with the primary structure.

Caveat worth knowing: PeSTo is partner-agnostic. Given a structure of the target
bound to something else, it will faithfully report *that* interface. Aggregating
across structures therefore aggregates across binding partners -- only pass in
structures whose other chains are irrelevant to the interaction you are docking.
"""
import argparse, os, sys


def ca_records(pdb):
    """CA lines in file order as (resnum, bfactor, line_index)."""
    out = []
    for i, line in enumerate(open(pdb)):
        if line.startswith(("ATOM", "HETATM")) and line[12:16].strip() == "CA":
            try:
                out.append((int(line[22:26]), float(line[60:66]), i))
            except ValueError:
                pass
    return out


def mapped_scores(pesto_pdb, src_pdb):
    """{resnum in src numbering: PeSTo score}, or None if the pair does not line up."""
    src = ca_records(src_pdb)
    out = ca_records(pesto_pdb)
    if not src or len(src) != len(out):
        return None
    return {s[0]: o[1] for s, o in zip(src, out)}


def aggregate(per, strategy, thr):
    allres = sorted({r for d in per.values() for r in d})
    mat = {r: [d[r] for d in per.values() if r in d] for r in allres}
    n = len(per)
    if strategy == "union":
        return {r for r in allres if max(mat[r]) >= thr}
    if strategy == "consensus":
        return {r for r in allres if sum(x >= thr for x in mat[r]) >= max(1, n // 2)}
    if strategy == "mean":
        return {r for r in allres if sum(mat[r]) / len(mat[r]) >= thr}
    if strategy == "best":
        k = max(per, key=lambda k: max(per[k].values()))
        return {r for r, x in per[k].items() if x >= thr}
    raise ValueError(f"unknown strategy {strategy}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--primary_pesto", required=True, help="PeSTo output for the docked structure")
    ap.add_argument("--primary_input", required=True, help="the structure that PeSTo output came from")
    ap.add_argument("--extra", nargs="*", default=[],
                    help="pairs of <pesto_output>:<its_input> for the additional structures")
    ap.add_argument("--strategy", default="union",
                    choices=["union", "consensus", "mean", "best"])
    ap.add_argument("--threshold", type=float, default=0.3)
    ap.add_argument("--out", required=True, help="PeSTo-format PDB to write")
    a = ap.parse_args()

    per = {}
    prim = mapped_scores(a.primary_pesto, a.primary_input)
    if prim is None:
        print(f"Error: {a.primary_pesto} and {a.primary_input} have different CA counts.")
        sys.exit(1)
    per["PRIMARY"] = prim
    for spec in a.extra:
        if ":" not in spec:
            print(f"  skip malformed --extra entry (want out:in): {spec}"); continue
        po, pi = spec.split(":", 1)
        if not (os.path.exists(po) and os.path.exists(pi)):
            print(f"  skip missing: {spec}"); continue
        sc = mapped_scores(po, pi)
        if sc is None:
            print(f"  skip {os.path.basename(pi)}: CA count differs from its PeSTo output")
            continue
        per[os.path.basename(pi)] = sc

    act = aggregate(per, a.strategy, a.threshold)
    prim_only = {r for r, x in prim.items() if x >= a.threshold}
    print(f"sources: {len(per)}  strategy: {a.strategy}  threshold: {a.threshold}")
    print(f"  primary alone : {len(prim_only)} active residues "
          f"(max score {max(prim.values()):.3f})")
    print(f"  aggregated    : {len(act)} active residues "
          f"({len(act - prim_only):+d} new, {len(prim_only - act)} dropped)")
    if a.strategy != "union" and (prim_only - act):
        print("  NOTE: this strategy discarded residues the primary structure had. "
              "'union' cannot do that.")

    # rewrite the primary PeSTo file, keeping its own numbering, with aggregated
    # residues set to 1.00 and everything else to 0.00
    src = ca_records(a.primary_input)
    lines = open(a.primary_pesto).readlines()
    k = 0
    for i, line in enumerate(lines):
        if line.startswith(("ATOM", "HETATM")) and line[12:16].strip() == "CA":
            resnum = src[k][0] if k < len(src) else None
            k += 1
            lines[i] = line[:60] + f"{(1.0 if resnum in act else 0.0):6.2f}" + line[66:]
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    open(a.out, "w").writelines(lines)
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
