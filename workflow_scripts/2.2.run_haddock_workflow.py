import subprocess
import os
import shutil
import argparse
import sys

# ==========================================
# Configuration Parameters
# ==========================================
# Fixed parameters for CRBN (Partner 2)
CRBN_PDB = "CRBN_pesto.pdb"
# ==========================================

def get_chain_id(pdb_file):
    """Reads the first chain ID from a PDB file."""
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM"):
                # Chain ID is at column 21 (0-indexed)
                return line[21]
    return None

def renumber_chain(input_pdb, output_pdb, new_chain_id):
    """Creates a new PDB file with the specified chain ID."""
    print(f"Renumbering {input_pdb} to chain {new_chain_id} -> {output_pdb}")
    with open(input_pdb, 'r') as fin, open(output_pdb, 'w') as fout:
        for line in fin:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Replace character at index 21 with new_chain_id
                # Careful with string immutability
                new_line = line[:21] + new_chain_id + line[22:]
                fout.write(new_line)
            else:
                fout.write(line)

def generate_config(pdb1_path, pdb2_path, restraint_path, run_dir, ncores, output_config, mode="mpi", iniseed=None):
    iniseed_line = f"iniseed = {iniseed}" if iniseed is not None else ""
    """Generates the HADDOCK3 configuration file."""
    config_content = f"""# ====================================================================
# Protein-protein docking example with NMR-derived ambiguous interaction restraints

# directory in which the scoring will be done
run_dir = "{run_dir}"

# execution mode
mode = "{mode}"
#  mpi: multi-node (2 nodes x 96 tasks = ncores = 192)
#  local: single node, multiprocessing with `ncores` workers
ncores = {ncores}

debug = true

# molecules to be docked
molecules =  [
    "{pdb1_path}",   
    "{pdb2_path}"
    ]

# ====================================================================
# Parameters for each stage are defined below, prefer full paths
# ====================================================================
[topoaa]
autohis = false

[rigidbody]
tolerance = 5
ambig_fname = "{restraint_path}"
sampling = 1000
{iniseed_line}

[caprieval]
reference_fname = "{pdb1_path}"

[seletop]
select = 200

[caprieval]
reference_fname = "{pdb1_path}"

[flexref]
tolerance = 5
ambig_fname = "{restraint_path}"

[caprieval]
reference_fname = "{pdb1_path}"

[mdref]
tolerance = 5
ambig_fname = "{restraint_path}"

[caprieval]
reference_fname = "{pdb1_path}"

[clustfcc]

[seletopclusts]
top_models = 4

[caprieval]
reference_fname = "{pdb1_path}"

# ====================================================================
"""
    with open(output_config, 'w') as f:
        f.write(config_content)
    print(f"Generated configuration file: {output_config}")

def main():
    parser = argparse.ArgumentParser(description="Run HADDOCK3 workflow.")
    parser.add_argument("--target_pdb", required=True, help="Path to the target PDB file (Partner 1).")
    parser.add_argument("--restraint", default="restraint.tbl", help="Path to the restraint file.")
    parser.add_argument("--run_dir", default="run1-mdref-full", help="HADDOCK3 run directory.")
    parser.add_argument("--config", default="docking.cfg", help="Output configuration filename.")
    parser.add_argument("--ncores", type=int, default=48, help="Number of cores to use.")
    parser.add_argument("--iniseed", type=int, default=None,
                        help="CNS random seed for the rigidbody stage. Omitted, "
                             "haddock3 uses its own default (917) and the run "
                             "behaves exactly as before. Note that fixing it "
                             "does not by itself make HADDOCK reproducible: "
                             "the parallel CNS workers finish in whatever order "
                             "they finish, and with sampling=1000 that reorders "
                             "the clustering.")
    parser.add_argument("--mode", default="mpi", choices=["mpi", "local"],
                        help="HADDOCK3 execution mode. 'local' uses multiprocessing on one "
                             "node and avoids the CUDA-enabled OpenMPI stack.")
    
    args = parser.parse_args()
    
    pdb1 = CRBN_PDB
    pdb2 = args.target_pdb
    
    # Check if CRBN_pesto.pdb exists
    if not os.path.exists(pdb1):
        print(f"Error: Fixed partner PDB {pdb1} not found in current directory.")
        sys.exit(1)

    # 1. Enforce Chain IDs: Target -> A, Partner -> B
    chain1 = get_chain_id(pdb1)
    chain2 = get_chain_id(pdb2)
    
    print(f"Original Chain ID for {pdb1}: {chain1}")
    print(f"Original Chain ID for {pdb2}: {chain2}")
    
    final_pdb1 = pdb1
    final_pdb2 = pdb2
    
    # Enforce pdb1 -> Chain A
    if chain1 != 'A':
        print(f"Renumbering {pdb1} to Chain A...")
        base_name1 = os.path.basename(pdb1)
        final_pdb1 = f"renumbered_A_{base_name1}"
        renumber_chain(pdb1, final_pdb1, 'A')
        
    # Enforce pdb2 -> Chain B
    if chain2 != 'B':
        print(f"Renumbering {pdb2} to Chain B...")
        base_name2 = os.path.basename(pdb2)
        final_pdb2 = f"renumbered_B_{base_name2}"
        renumber_chain(pdb2, final_pdb2, 'B')
    
    # 2. Generate Configuration
    generate_config(final_pdb1, final_pdb2, args.restraint, args.run_dir,
                    args.ncores, args.config, args.mode, args.iniseed)
    
    # 3. Run HADDOCK3
    if os.path.exists(args.run_dir):
        print(f"Run directory {args.run_dir} exists. Running with --restart 0.")
        command = ["haddock3", args.config, "--restart", "0"]
    else:
        print(f"Run directory {args.run_dir} does not exist. Running normally.")
        command = ["haddock3", args.config]

    print(f"Running command: {' '.join(command)}")
    
    try:
        subprocess.run(command, check=True)
        print("HADDOCK3 run completed successfully.")
    except subprocess.CalledProcessError as e:
        print("HADDOCK3 run failed.")
        print(e)
        sys.exit(1)

if __name__ == "__main__":
    main()
