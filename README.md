# GluePlex

GluePlex is an automated pipeline for molecular glue docking and analysis, integrating **PeSTo** for interface prediction, **HADDOCK3** for docking, and **Boltz2** for confidence scoring and refinement.

The `GluePlex_doublepass` branch adds a **Step 4** that re-runs Boltz2 using the top predictions from Step 3 as multi-chain templates, concentrating sampling around the best interfaces, then produces a combined analysis of both passes.

## Overview

The pipeline consists of four steps:
1.  **Step 1: PeSTo** — Predicts active residues and interfaces on the target protein.
2.  **Step 2: HADDOCK3** — Performs data-driven docking using the predicted interfaces to generate initial models.
3.  **Step 3: Boltz2 (Pass 1)** — Refines the HADDOCK3 models and scores them using ipTM confidence metrics, followed by PCA/UMAP visualization.
4.  **Step 4: Boltz2 (Pass 2)** — Selects the top-N models from Pass 1 (ranked by `pair_chains_ipTM_avg`), converts them to mmCIF templates, re-runs Boltz2 with multi-chain templating, then merges both passes for unified post-analysis and visualization.

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

## Configuration

The pipeline is controlled by `project_config.yaml`. Edit this file to define your inputs and parameters.

```yaml
# General Settings
base_dir: "/path/to/GluePlex"

# Step 1: PeSTo
target_pdb: "./VAV1_af2.pdb"  # Path to target PDB or PDB ID (auto parsing e.g., 1E9H_A)

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
