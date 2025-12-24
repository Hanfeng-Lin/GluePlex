import subprocess
import os
import shutil

# Configuration
PDB1 = "CRBN_pesto.pdb"
PDB2 = "6S76_pesto.pdb"
RESTRAINT_FILE = "restraint.tbl"
RUN_DIR = "run1-mdref-full"
CONFIG_FILE = "docking-protein-protein-mdref-mpi-full.cfg"
NCORES = 48

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

def generate_config(pdb1_path, pdb2_path, restraint_path, output_config):
    """Generates the HADDOCK3 configuration file."""
    config_content = f"""# ====================================================================
# Protein-protein docking example with NMR-derived ambiguous interaction restraints

# directory in which the scoring will be done
run_dir = "{RUN_DIR}"

# execution mode 
mode = "mpi"
#  2 nodes x 96 tasks = ncores = 192
ncores = {NCORES}

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
autohis = true

[rigidbody]
tolerance = 5
ambig_fname = "{restraint_path}"
sampling = 1000

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
    # 1. Check Chain IDs
    chain1 = get_chain_id(PDB1)
    chain2 = get_chain_id(PDB2)
    
    print(f"Chain ID for {PDB1}: {chain1}")
    print(f"Chain ID for {PDB2}: {chain2}")
    
    final_pdb2 = PDB2
    
    if chain1 == chain2:
        print("Chain IDs conflict! Renumbering PDB2...")
        # Assume we want to change PDB2 to chain 'B' if it's 'A', or 'A' if it's 'B'
        new_chain = 'B' if chain2 == 'A' else 'A'
        final_pdb2 = f"renumbered_{PDB2}"
        renumber_chain(PDB2, final_pdb2, new_chain)
    
    # 2. Generate Configuration
    generate_config(PDB1, final_pdb2, RESTRAINT_FILE, CONFIG_FILE)
    
    # 3. Run HADDOCK3
    if os.path.exists(RUN_DIR):
        print(f"Run directory {RUN_DIR} exists. Running with --restart 0.")
        command = ["haddock3", CONFIG_FILE, "--restart", "0"]
    else:
        print(f"Run directory {RUN_DIR} does not exist. Running normally.")
        command = ["haddock3", CONFIG_FILE]

    print(f"Running command: {' '.join(command)}")
    
    try:
        subprocess.run(command, check=True)
        print("HADDOCK3 run completed successfully.")
    except subprocess.CalledProcessError as e:
        print("HADDOCK3 run failed.")
        print(e)

if __name__ == "__main__":
    main()
