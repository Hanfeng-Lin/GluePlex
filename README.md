# GluePlex

GluePlex is an automated pipeline for molecular glue docking and analysis, integrating **PeSTo** for interface prediction, **HADDOCK3** for docking, and **Boltz2** for confidence scoring and refinement.

**Step 4** re-runs Boltz2 using the top predictions from Step 3 as multi-chain templates, concentrating sampling around the best interfaces, then produces a combined analysis of both passes.

The pipeline can be run in six **arms** that differ in how the docking template is built, from the full pipeline down to Boltz2 on its own. See [Choosing an arm](#choosing-an-arm).

## Overview

The pipeline consists of four steps, each orchestrated by `run_glueplex.py` and implemented by numbered scripts in `workflow_scripts/`:

1.  **Step 1: PeSTo** — Predicts interface residues on the target protein.
2.  **Step 2: HADDOCK3** — Data-driven docking using the predicted interface to generate initial ternary models.
3.  **Step 3: Boltz2 (Pass 1)** — Re-predicts each HADDOCK cluster with Boltz2, scores by ipTM, and visualizes with PCA/UMAP.
4.  **Step 4: Boltz2 (Pass 2)** — Re-runs Boltz2 using the top-N Pass-1 models as multi-chain templates, then merges both passes for unified analysis.

The sections below detail what each step actually does internally.

### Step 1 — PeSTo interface prediction

*Script: `tools/PeSTo/apply_model.py` (run via the orchestrator)*

- Accepts `target_pdb` as either a local PDB file or a PDB ID (e.g. `1E9H_A`, with optional chain selection). Files are copied into `input/`; IDs are fetched on the fly.
- Runs the **PeSTo** deep-learning model (PyTorch) to score every residue for the probability of being part of a protein–protein interface.
- Writes the per-residue interface probability into the **B-factor column** of `step1_pesto/{base}_pesto.pdb`. This annotated PDB feeds Step 2.

### Step 2 — HADDOCK3 docking

Three sub-steps, all run in the `haddock3` environment:

- **2.0 `fetch_active_res.py`** — Derives the target's active (interface) residues when `target_active_res` is not set in the config. It runs **NACCESS** to compute relative solvent accessibility (keeping residues >15% exposed), parses the PeSTo B-factors (keeping CA atoms with score ≥ 0.3), and takes the **intersection** (interface ∧ surface-exposed). It falls back to PeSTo-only or SASA-only if the intersection is empty, and prints the residues as a comma-separated list.
- **2.1 `generate_actpass.py`** — Builds HADDOCK3 ambiguous interaction restraints (`restraint.tbl`). The fixed partner (CRBN) uses hardcoded active residues `105,332,341`; the target uses the residues from 2.0. For each partner it derives passive residues via `haddock3-restraints passive_from_active`, then combines both into ambiguous restraints with `active_passive_to_ambig`.
- **2.2 `run_haddock_workflow.py`** — Normalizes chain IDs (fixed partner → chain A, target → chain B), generates a HADDOCK3 config, and runs `haddock3` in MPI mode (`ncores` from config). The docking protocol is: `topoaa` → `rigidbody` (sampling 1000) → `seletop` (top 200) → `flexref` → `mdref` → `clustfcc` → `seletopclusts` (top 4 models/cluster), with `caprieval` evaluation after each refinement stage.
- **2.3 `export_mmcif.py`** — Locates the `*_seletopclusts` directory, takes the top model of each cluster (`cluster_*model_1.pdb.gz`), and converts each to **mmCIF** with ChimeraX (`bestGuess true`) so Boltz2 can use them as templates.

### Step 3 — Boltz2 Pass 1

Two sub-steps, run in the `boltz_analysis` environment:

- **3.1 `run_batch_clusters.py`** — For each cluster mmCIF from Step 2, generates a Boltz2 input YAML (chain A = protein A sequence, chain B = protein B sequence, ligand C = SMILES, `msa: empty`, and a **multi-chain template** spanning chains A+B) plus a wrapper run script, then runs `boltz predict` with `diffusion_samples` samples. CIFs are distributed round-robin across the configured GPUs via a thread pool (one prediction per GPU at a time).
- **3.2 `post_analysis.py`** — Flattens all predicted PDBs/JSONs into `{protein_a}-{protein_b}/`, renumbers residues to match the configured ranges, then **superposes every model on chain A** (gemmi) and extracts chain B CA coordinates. It reads `ipTM` and `pair_chains_iptm` from the Boltz confidence JSONs into `atom.csv`, optionally computes reference RMSD (see `--ref`), writes `top10_models.log`, runs **PCA and UMAP** visualizations (interactive Bokeh HTML + static PNG), and packages the top-10 PDBs + plots into a `step3_top10_output/` folder.

### Step 4 — Boltz2 Pass 2 (double-pass refinement)

Three sub-steps, run in the `boltz_analysis` environment:

- **4.1 `prepare_second_pass.py`** — Reads Step 3's `atom.csv`, sorts the Pass-1 rows by `rank_metric` (default `pair_chains_ipTM_avg`), and selects the top-N. Each selected PDB is converted to mmCIF via ChimeraX (Boltz's PDB template parser fails on bare-ATOM PDBs without SEQRES), and a per-template Boltz YAML + wrapper script is written into `step4_boltz/`. Emits `pass2_manifest.txt` and `pass2_template_mapping.log` (source → renamed template). Does not run Boltz itself.
- **4.2 `run_second_pass.py`** — Reads the manifest and launches every pass-2 Boltz job in parallel, with concurrency bounded by the number of unique GPUs the jobs target.
- **4.3 `merge_and_analyze.py`** — Symlinks all Pass-1 and Pass-2 PDBs/JSONs into `merged_pdbs/`, then re-runs `3.2.post_analysis.py` on the merged set to produce `atom_2ndpass.csv` and combined PCA/UMAP plots (Pass-2 points are visually distinguished — see [Analysis](#analysis)).

## Prerequisites

The pipeline relies on three specific Conda environments being available on your system. Please ensure these are set up before running the orchestrator:

1.  **`pesto`**: Environment containing the PeSTo model and dependencies (PyTorch, etc.).
2.  **`haddock3`**: Environment containing the HADDOCK3 executable.
3.  **`boltz_analysis`**: Environment containing Boltz2 and analysis tools (Gemmi, BioPython, scikit-learn, Bokeh, etc.).

*Note: The orchestrator (`run_glueplex.py`) explicitly calls `conda run -n [env_name]`. Ensure your environments match these names or update the script/config accordingly.*

## Installation

1.  **Clone the repository**:
    ```bash
    git clone https://github.com/Hanfeng-Lin/GluePlex.git
    cd GluePlex
    git checkout GluePlex_doublepass
    ```

2.  **Install NACCESS and UCSF ChimeraX**:
    Ensure `naccess` is installed and available in your PATH. This is required for calculating surface accessibility in Step 2.
    ChimeraX is required for converting PDB to mmCIF in Step 2 and Step 4.

    ```bash
    which naccess
    which chimerax
    ```
    https://www.cgl.ucsf.edu/chimerax/download.html

3.  **Set up PeSTo**:
    Since PeSTo cannot be distributed directly, you need to clone it into the `tools` directory:
    ```bash
    cd tools
    git clone https://github.com/LBM-EPFL/PeSTo.git
    cp apply_model.py PeSTo/
    cd ..
    ```

4.  **Create the required Conda environments**:
    
    **PeSTo**
    ```bash
    conda env create -f env_pesto.yml
    ```

    **HADDOCK3**
    ```bash
    conda env create -f env_haddock3.yml
    ```

    **Boltz Analysis**
    ```bash
    conda env create -f env_boltz_analysis.yml
    ```

5.  **Orchestrator Environment** (Optional)
    The orchestrator (`run_glueplex.py`) requires `PyYAML`. You can install it in your base environment or activate the pesto environment:
    ```bash
    conda activate pesto
    ```

## Choosing an arm

The same pipeline can be run six ways. They differ in exactly one thing — how
much structural guidance Boltz2 is given before it predicts — so running more
than one tells you where a result is coming from.

Pick one by copying the matching file out of `examples/`. Five of the six are
configuration only and use the normal entry point; the sixth has its own
runner because it skips docking entirely.

| arm | what it docks | interface residues | config | command |
|---|---|---|---|---|
| `af` | an AlphaFold model | PeSTo on that model | `examples/arm_af.yaml` | `run_glueplex.py` |
| `aggaf` | an AlphaFold model | PeSTo unioned over **every** structure you supply | `examples/arm_aggaf.yaml` | `run_glueplex.py` |
| `medoid` | the most central experimental chain | same union | `examples/arm_medoid.yaml` | `run_glueplex.py` |
| `pesto` | the experimental chain PeSTo scores highest | same union | `examples/arm_pesto.yaml` | `run_glueplex.py` |
| `nopesto` | an AlphaFold model | **none** — HADDOCK falls back to SASA | `examples/arm_nopesto.yaml` | `run_glueplex.py` |
| `msa` | nothing — no docking at all | — | `examples/arm_msa.yaml` | `run_msa_baseline.py` |

```bash
cp examples/arm_af.yaml my_target.yaml
# edit: target_pdb, the two sequences, ligand_smiles, gpus
python run_glueplex.py --config my_target.yaml
```

### Which one do you want?

**Start with `af`.** It needs one AlphaFold model of your target and nothing
else, and it is the arm the others are measured against.

**Use `aggaf`, `medoid` or `pesto` if your target has solved structures.** All
three run PeSTo on every structure you list in `pesto_extra_pdbs` and union the
predicted interface residues; they differ only in which coordinates get docked.
`aggaf` still docks the AlphaFold model, so comparing it against `af` isolates
what the extra interface evidence is worth, and comparing `medoid` or `pesto`
against `aggaf` isolates what the experimental coordinates are worth. Compare
`medoid` or `pesto` against `af` directly and you cannot tell those two apart,
because both changed at once.

`medoid` picks the chain with the lowest mean CA RMSD to all the others, which
needs no scoring model. `pesto` picks the chain PeSTo is most confident about.
Neither is reliably better; they disagree often enough to be worth running both
if the target matters.

**Use `nopesto` to find out whether the interface prediction is doing anything.**
Step 1 is skipped and the target is written out with its B-factor column zeroed,
so `2.0.fetch_active_res.py` finds no PeSTo residues and takes its SASA branch —
the docking is restrained to the whole solvent-accessible surface. CRBN is
untouched and keeps its hand-picked actives, so this removes the *target's*
interface prediction, not the restraints as a whole. Needs NACCESS installed.

**Use `msa` as the floor.** Off-the-shelf Boltz2: two sequences, a SMILES, an
MSA from the MMseqs2 server, no PeSTo, no HADDOCK, no template. This is what
published cofolding benchmarks measure, and the gap between it and any other
arm is what the docking machinery is worth. It needs no structure of your
target at all, so it is also the only arm you can run on a target that has none.

```bash
conda activate boltz_analysis          # `boltz` must be on PATH
cp examples/arm_msa.yaml my_target.yaml
python run_msa_baseline.py --config my_target.yaml --samples 440
```

`--samples` defaults to 440 to match the median number of models a full
two-pass run produces, so the comparison is at equal sampling budget rather
than at equal nominal settings. Sampling is chunked because Boltz draws a whole
request in one batch and a large one runs out of memory. `--msa-cache DIR`
reuses server MSAs across runs sharing a sequence.

### Keeping arms comparable

`pesto_aggregate` must stay at `union` in every arm that aggregates. It is the
constant that makes the three aggregating arms differ only in which structure
is docked; changing it introduces a second variable and the arms can no longer
be read against each other. The other three merge strategies the aggregator
accepts (`consensus`, `mean`, `best`) are there for a separate experiment, not
for mixing into this comparison.

Run each arm into its own `project_name`, or the later one will reuse the
earlier one's step directories.


## Configuration

The pipeline is controlled by a YAML config passed with `--config`. Start from
the arm you want in `examples/` rather than from `project_config.yaml`, which is
wired to one specific run on one host. The keys below are the full set; an arm's
example file is this with the irrelevant parts removed.

```yaml
# General Settings
base_dir: "/path/to/GluePlex"

# Step 1: PeSTo
target_pdb: "./VAV1_af2.pdb"  # Path to target PDB or PDB ID (auto parsing e.g., 1E9H_A)

# Optional. Other structures of the same target; PeSTo runs on each and the
# predicted interface residues are merged into target_pdb's numbering.
# Unset = single-structure behaviour, i.e. the `af` arm.
# pesto_extra_pdbs:
#   - "./structures/TARGET_1abc_A.pdb"
# pesto_aggregate: union        # union | consensus | mean | best — keep at union
# pesto_threshold: 0.3          # a residue is active at or above this score

# Optional. Skip PeSTo and write target_pdb out with zeroed B-factors, so
# step 2 finds no interface residues and falls back to SASA. The `nopesto`
# arm. Requires target_pdb to be a local file, not an ID.
# skip_pesto: true

# Step 2: HADDOCK3
# Active residues for Target. If empty, they will be fetched automatically
# through PeSTo score and SASA filtering.
# target_active_res: "1,3,9,10,..."
haddock_run_dir: "haddock_results"
haddock_ncores: 80

# Step 3: Boltz2 (Pass 1)
step3:
  protein_a_name: "CRBN"
  protein_a_sequence: "..."
  protein_b_name: "VAV1"
  protein_b_sequence: "..."
  ligand_name: "NGT-201-12"
  ligand_smiles: "..."
  diffusion_samples: 20
  gpus: [0,1,2,3,4,5,6,7]

# Step 4: Boltz2 (Pass 2 — double-pass refinement)
step4:
  top_n: 10                        # number of pass-1 templates to carry into pass-2
  diffusion_samples: 20            # samples per pass-2 job
  rank_metric: pair_chains_ipTM_avg
  gpus: [0,1,2,3,4,5,6,7]
```

## Usage

Run the pipeline using the `run_glueplex.py` script.

### Run the Full Pipeline (all 4 steps)
```bash
python run_glueplex.py --config project_config.yaml
```

The project directory is automatically named `Result_{protein_b}_{protein_a}_GluePlex2pass` (e.g., `Result_VAV1_CRBN_GluePlex2pass`).

### Run without double-pass (Steps 1-3 only)
```bash
python run_glueplex.py --config project_config.yaml --no-doublepass
```

### Run Specific Steps
You can run individual steps using the `--step` argument:

*   **Step 1 only (PeSTo)**:
    ```bash
    python run_glueplex.py --step 1
    ```
*   **Step 2 only (HADDOCK3)**:
    ```bash
    python run_glueplex.py --step 2
    ```
*   **Step 3 only (Boltz2 Pass 1)**:
    ```bash
    python run_glueplex.py --step 3
    ```
*   **Step 4 only (Boltz2 Pass 2)**:
    ```bash
    python run_glueplex.py --step 4
    ```

### Override Existing Results
To force re-running Boltz predictions even if outputs exist:
```bash
python run_glueplex.py --step 3 --override
```

### Reference PDB for RMSD Evaluation

When a ground-truth crystal structure is available, use `--ref` to compute chain B CA RMSD for each predicted model:

```bash
python run_glueplex.py --config project_config.yaml --step 3 --ref groundtruth.pdb
```

The reference PDB must contain exactly 3 chains (A: fixed partner, B: target, C: ligand). For each model, chain A is aligned to the reference chain A using sequence alignment (BioPython `PairwiseAligner`) and structural superposition (`Superimposer`), then chain B CA RMSD is measured against the reference chain B. This handles different residue numbering, chain truncations, and missing density gaps.

The `ref_rmsd` column is added to `atom.csv` / `atom_2ndpass.csv`, and an additional `PCA_ref_rmsd` / `UMAP_ref_rmsd` plot is generated (color range 0–10 Å, clamped).

For targets with multiple ground-truth conformations (e.g., active/inactive), the helper script can be run directly:
```bash
# Compute RMSD against active conformation
python workflow_scripts/_ref_align_rmsd.py \
    --ref_pdb active.pdb \
    --input_csv atom.csv \
    --results_dir CRBN-TARGET \
    --align_chain A --extract_chain B

# Regenerate plots, then rename with suffix
python workflow_scripts/_pca_coord_energy.py \
    --input_csv atom.csv \
    --ligand_name LIGAND --protein_name TARGET --protein_a_name CRBN
mv PCA_ref_rmsd_CRBN_TARGET.png PCA_ref_active_CRBN_TARGET.png
```

## Output Structure

```
Result_VAV1_CRBN_GluePlex2pass/
├── input/                          # Copies of input files and config
├── logs/                           # Log files for each step
│   ├── step1_pesto.log
│   ├── step2_haddock.log
│   ├── step3_boltz.log
│   └── step4_double_pass.log
├── step1_pesto/                    # PeSTo outputs (predicted PDBs)
├── step2_haddock/                  # HADDOCK3 run directory and mmCIF exports
│   ├── haddock_results/            # Full HADDOCK3 output
│   ├── restraint.tbl               # Generated restraints
│   └── *.cif                       # Exported models for Boltz
├── step3_boltz/                    # Boltz Pass 1 predictions and analysis
│   ├── boltz_results/              # Raw Boltz outputs
│   ├── CRBN-VAV1/                  # Flattened results (PDBs, JSONs)
│   ├── atom.csv                    # Summary CSV with scores, coords, and ref_rmsd
│   ├── top10_models.log            # Top 10 models (with ref_rmsd if --ref used)
│   ├── PCA_*_CRBN_VAV1.html/png    # Interactive + static PCA plots
│   └── UMAP_*_CRBN_VAV1.html/png   # Interactive + static UMAP plots
└── step4_boltz/                    # Boltz Pass 2 (double-pass refinement)
    ├── templates_cif/              # Pass-1 PDBs converted to mmCIF for templating
    ├── pass2_manifest.txt          # List of pass-2 run scripts
    ├── pass2_template_mapping.log  # Source PDB → pass-2 template name mapping
    ├── boltz_results/              # Raw Boltz Pass 2 outputs
    ├── merged_pdbs/                # Symlinks to all Pass 1 + Pass 2 PDBs/JSONs
    ├── CRBN-VAV1_2ndpass/          # Flattened merged results
    ├── atom_2ndpass.csv            # Combined CSV for both passes
    ├── PCA_*_CRBN_VAV1.html/png    # Combined PCA plots
    └── UMAP_*_CRBN_VAV1.html/png   # Combined UMAP plots
```

## Analysis

### Plot Types

The post-analysis generates plot variants in both PCA and UMAP, each as interactive HTML (Bokeh) and static PNG:

- **`cluster`** — colors by HADDOCK3 cluster ID (Pass 1) or pass-2 template rank.
- **`chainpair_ipTM`** — colors by the average ipTM score between chain A (e.g., CRBN) and chain B (target protein), the most relevant metric for molecular glue interfaces.
- **`ipTM`** — colors by the overall ipTM score of the model.
- **`ref_rmsd`** — colors by chain B CA RMSD against the reference crystal structure (range 0–10 Å). Only generated when `--ref` is provided.

### Pass Distinction

In combined Pass 1 + Pass 2 plots, **Pass-2 points are drawn with black edges** to visually distinguish them from Pass-1 points, which have no edges. The `pass` column in `atom_2ndpass.csv` (1 or 2) encodes this programmatically.

### Reference RMSD

The `--ref` flag enables ground-truth comparison. The alignment uses:
1. **Sequence alignment** (`Bio.Align.PairwiseAligner`) to pair residues regardless of numbering
2. **Structural superposition** (`Bio.PDB.Superimposer`) on paired chain A CAs
3. **Chain B RMSD** measured after chain A alignment (no re-superposition of chain B)

This handles common mismatches between crystal structures and predicted models (different residue ranges, missing density, truncated chains).

### Viewing Plots on a Remote Server

```bash
cd Result_VAV1_CRBN_GluePlex2pass/step4_boltz
python -m http.server 8000
```
Then access `http://localhost:8000` in your browser (requires port forwarding).

### Key Log Files

- `step3_boltz/top10_models.log` — ranked list of the best Pass-1 predictions (with ref_rmsd if `--ref` was used).
- `step4_boltz/pass2_template_mapping.log` — shows which Pass-1 model became which Pass-2 template (source → rename mapping).

## Why Double-Pass?

Boltz2 uniquely supports **multi-chain templates** — a single template can encode the relative orientation of two protein chains (e.g., CRBN and the target). AlphaFold 3, by contrast, only supports single-chain templates and rejects multi-chain mmCIFs. By using top-scoring Pass-1 ternary complexes as templates, Pass 2 concentrates Boltz sampling around the most plausible CRBN–glue–target geometries, often improving interface confidence scores.
