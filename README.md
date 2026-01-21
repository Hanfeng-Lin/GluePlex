# GluePlex

GluePlex is an automated pipeline for molecular glue docking and analysis, integrating **PeSTo** for interface prediction, **HADDOCK3** for docking, and **Boltz** for confidence scoring and refinement.

## Overview

The pipeline consists of three main steps:
1.  **Step 1: PeSTo**: Predicts active residues and interfaces on the target protein.
2.  **Step 2: HADDOCK3**: Performs data-driven docking using the predicted interfaces to generate initial models.
3.  **Step 3: Boltz**: Refines the HADDOCK3 models and scores them using AlphaFold3-like confidence metrics (ipTM), followed by clustering and visualization.

## Prerequisites

The pipeline relies on three specific Conda environments being available on your system. Please ensure these are set up before running the orchestrator:

1.  **`pesto`**: Environment containing the PeSTo model and dependencies (PyTorch, etc.).
2.  **`haddock3`**: Environment containing the HADDOCK3 executable.
3.  **`boltz_analysis`**: Environment containing Boltz and analysis tools (Gemmi, UMAP, etc.).

*Note: The orchestrator (`run_glueplex.py`) explicitly calls `conda run -n [env_name]`. Ensure your environments match these names or update the script/config accordingly.*

## Installation

1.  **Clone the repository**:
    ```bash
    git clone <repository_url>
    cd GluePlex
    ```

2.  **Create the required Conda environments**:
    
    **Step 1: PeSTo**
    ```bash
    conda env create -f env_pesto.yml
    ```

    **Step 2: HADDOCK3**
    ```bash
    conda env create -f env_haddock3.yml
    ```
    *Note: You also need `chimerax` in your PATH for mmCIF conversion.*

    ```bash
    which chimerax
    ```

    **Step 3: Boltz Analysis**
    ```bash
    conda env create -f env_boltz_analysis.yml
    ```

3.  **Orchestrator Environment** (Optional)
    The orchestrator (`run_glueplex.py`) requires `PyYAML`. You can install it in your base environment or enable pesto environment:
    ```bash
    conda activate pesto
    ```

## Configuration

The pipeline is controlled by `project_config.yaml`. Edit this file to define your inputs and parameters.

```yaml
# General Settings
project_name: "Result_VAV1_CRBN"
base_dir: "/path/to/GluePlex"

# Step 1: PeSTo
target_pdb: "./VAV1_af2.pdb"  # Path to target PDB or PDB ID (auto parsing e.g., 1E9H_A)

# Step 2: HADDOCK3
# Active residues for Target (Partner 1). 
# If empty, they will be fetched automatically (TODO: SASA filtering).
target_active_res: "1,3,9,10,..." 
haddock_run_dir: "haddock_results"
haddock_ncores: 80

# Step 3: Boltz
step3:
  protein_a_name: "CRBN"
  protein_a_sequence: "..."
  protein_b_name: "VAV1"
  protein_b_sequence: "..."
  ligand_name: "NGT-201-12"
  ligand_smiles: "..."
  diffusion_samples: 20
  gpus: [0,1,2,3,4,5,6,7] # List of GPU IDs to use
```

## Usage

Run the pipeline using the `run_glueplex.py` script.

### Run the Full Pipeline
```bash
python run_glueplex.py --config project_config.yaml
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
*   **Step 3 only (Boltz)**:
    ```bash
    python run_glueplex.py --step 3
    ```

### Override Existing Results
To force re-running Boltz predictions even if outputs exist:
```bash
python run_glueplex.py --step 3 --override
```

## Output Structure

The pipeline creates a project directory (defined in `project_config.yaml`) with the following structure:

```
Result_VAV1_CRBN/
├── input/                  # Copies of input files and config
├── logs/                   # Log files for each step
├── step1_pesto/            # PeSTo outputs (predicted PDBs)
├── step2_haddock/          # HADDOCK3 run directory and mmCIF exports
│   ├── haddock_results/    # Full HADDOCK3 output
│   ├── restraint.tbl       # Generated restraints
│   └── *.cif               # Exported models for Boltz
└── step3_boltz/            # Boltz predictions and analysis
    ├── boltz_results/      # Raw Boltz outputs
    ├── CRBN-VAV1/          # Flattened results (PDBs, JSONs)
    ├── atom.csv            # Summary CSV of scores and coordinates
    ├── UMAP_chainpair_ipTM.html # Interactive UMAP plot
    └── ...
```

## Analysis

The final step generates UMAP plots in png format and interactive plots (`.html` files) in the `step3_boltz` directory. These plots visualize the conformational space of the docked models, colored by their confidence scores (ipTM).

UMAP_cluster colors by cluster id from HADDOCK3.
UMAP_chainpair_ipTM colors by the average ipTM score between chain 0 (e.g. CRBN) and chain 1 (target protein, e.g. VAV1).
UMAP_ipTM colors by the overall ipTM score of the model.

To view the plots on a remote server, you can start a simple HTTP server:
```bash
cd Result_VAV1_CRBN/step3_boltz
python -m http.server 8000
```
Then access `http://localhost:8000` in your browser (requires port forwarding).
