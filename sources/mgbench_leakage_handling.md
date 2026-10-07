# MGBench / MG-PDB — how they handled seen cases

Verified against the published PDF (JCIM 2025, 65, 11136-11148,
doi 10.1021/acs.jcim.5c01860). Supersedes an earlier version of this file that
was written from secondhand summaries and got the homology criterion wrong.

## Three splits, not one

1. RESOLUTION. 221 noncovalent MG ternary complexes -> 202. XRD worse than
   3.5 A removed; EM structures kept if the LOCAL resolution at the MG-binding
   interface meets 3.5 A even when the overall value does not.

2. TEMPORAL. Release date after 2021-09-30 -> MGBench, 88 structures. The
   remaining 114 are the "seen set" and are REPORTED AS A COMPARISON GROUP in
   Figure 3, not silently dropped.

3. HOMOLOGY, within the 88. MMseqs2 against all PDB sequences released before
   2021-09-30, following AF3's evaluation-set interface criteria. Verbatim:

     "If both polymers have a length of at least 16 amino acids and exhibit
      greater than 40% sequence identity to two chains in the same complex in
      the training data, the structure is filtered out."
     "For structures containing chains with fewer than 16 residues, the
      similarity of the longer chain must be less than 40% sequence identity to
      any structure in the training data."

   -> low-homology 25, high-homology 63.

   NOTE the criterion is a PAIR test -- both chains matching the SAME training
   complex -- not each chain independently. A novel neosubstrate against a
   well-covered E3 stays in low-homology.

4. QUALITY (a fourth split, often overlooked). HQ 23 / LQ 65. The paper finds
   the HQ/LQ gap large, so quality confounds any homology comparison that does
   not control for it.

## Headline numbers (AF3, best method, on MGBench as a whole)

- PPI interface prediction 50.6% success
- MG-protein interaction recovery 32.9% (all methods ~30%)
- low-homology significantly below high-homology; the paper states the gap
  qualitatively and in Figure 3 rather than as a single number in the text.

## The limitation the paper itself documents — and it is about CRBN

Case study: CRBN:dWIZ-1:WIZ(ZF7) (PDB 8TZX). WIZ is in the LOW-homology set
(<40% identity) and yet every method modelled it well. The paper's explanation:
the C2H2 G-loop degron motif is structurally conserved, and CRBN/C2H2 ternary
structures ARE in the training set, naming

  5FQD CK1a-lenalidomide, 5HXB GSPT1-CC885, 6H0G ZNF692-pomalidomide,
  6UML SALL4-thalidomide, 7LPS IKZF2-ALV1

So sequence identity UNDERSTATES leakage for CRBN zinc-finger degrons: a target
can be sequence-novel and still be effectively seen through a conserved motif.
The paper also notes that nearly all domain-motif structures in MGBench are
high-homology, and that only CRBN and 14-3-3 systems were present in training.

## Not done by MGBench

- No assessment of LIGAND novelty; only protein sequence homology.
- No ligand truncation experiments.

## Relation to the GluePlex stratification

All five structures the paper names as being in the training set
(CSNK1A1, GSPT1, IKZF2, SALL4, ZNF692) fall in GluePlex's SEEN (direct)
stratum -- independent confirmation of that assignment.

WIZ is SEEN (homolog) in GluePlex but LOW-homology in MGBench. GluePlex is the
more conservative call on exactly the case the paper flags as a failure of the
40% rule.

GluePlex's novel stratum -- ALK, CDK2, ELAVL1, G3BP2, HBS1L, KAT2A, MTOR, NEK7,
PDE6D, PPIL4, VAV1 -- contains no C2H2 zinc finger, so it is motif-disjoint
from the G-loop class and not exposed to this failure mode.

Difference in scope: MGBench excludes pre-cutoff structures from the benchmark
and reports them separately as the "seen set"; GluePlex keeps all 81 molecules
and reports SEEN (direct) 34 / SEEN (homolog) 12 / novel 35 as strata of one
set. MGBench's "high-homology" is a homolog-level match within post-cutoff
structures, so it is closer to SEEN (homolog) than to SEEN (direct).
