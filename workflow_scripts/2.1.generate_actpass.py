import subprocess
import os
import argparse
import sys

# ==========================================
# Configuration Parameters
# ==========================================
# Fixed parameters for CRBN (Partner 1)
CRBN_PDB = "CRBN_pesto.pdb"
# CRBN_pesto.pdb numbering (1-397) == UniProt 46-442, i.e. UniProt = local + 45.
# Default kept for backward compatibility: F150, E377, W386.
# Validation against 26 crystal ternaries: W386 contacts the glue in 26/26 and the
# neosubstrate in 18/26; E377 22/26 and 10/26; but F150 contacts the glue in only
# 3/26 and sits in the N-terminal Lon domain. Override with --crbn_res, e.g. the
# glue-binding pocket 306,307,308,333,334,335,341,355,357 (UniProt 351,352,353,
# 378,379,380,386,400,402), which is derivable from binary CRBN-thalidomide
# structures and so leaks nothing about the neosubstrate.
CRBN_ACTIVE_RES = "105,332,341"
# ==========================================

def run_haddock3_restraints(pdb_file, active_residues, output_filename):
    """
    Runs haddock3-restraints passive_from_active and writes the output to a file.
    """
    command = ["haddock3-restraints", "passive_from_active", pdb_file, active_residues]
    
    print(f"Running command: {' '.join(command)}")
    
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        passive_residues = result.stdout.strip()
        
        # Format active residues: replace commas with spaces
        formatted_active = active_residues.replace(",", " ")
        
        with open(output_filename, "w") as f:
            f.write(f"{formatted_active}\n")
            f.write(f"{passive_residues}\n")
            
        print(f"Successfully generated {output_filename}")
        
    except subprocess.CalledProcessError as e:
        print(f"Error running command for {pdb_file}:")
        print(f"STDOUT: {e.stdout}")
        print(f"STDERR: {e.stderr}")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)

def run_haddock3_ambig(actpass_file1, actpass_file2, output_filename):
    """
    Runs haddock3-restraints active_passive_to_ambig and writes the output to a file.
    """
    # Assuming haddock3-restraints is in the path or using the full path if needed
    # Using just "haddock3-restraints" assuming the env is active
    command = ["haddock3-restraints", "active_passive_to_ambig", actpass_file1, actpass_file2]
    
    print(f"Running command: {' '.join(command)}")
    
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        restraints = result.stdout
        
        with open(output_filename, "w") as f:
            f.write(restraints)
            
        print(f"Successfully generated {output_filename}")
        
    except subprocess.CalledProcessError as e:
        print(f"Error running command for {output_filename}:")
        print(f"STDOUT: {e.stdout}")
        print(f"STDERR: {e.stderr}")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Generate HADDOCK3 restraints.")
    parser.add_argument("--target_pdb", required=True, help="Path to the target PDB file (Partner 1).")
    parser.add_argument("--target_res", required=True, help="Comma-separated active residues for the target.")
    parser.add_argument("--output", default="restraint.tbl", help="Output filename for the restraints table.")
    parser.add_argument("--crbn_res", default=None,
                        help="Comma-separated CRBN active residues in CRBN_pesto.pdb numbering. "
                             f"Default: {CRBN_ACTIVE_RES}")

    args = parser.parse_args()
    crbn_active = args.crbn_res or CRBN_ACTIVE_RES

    # 1. Generate actpass for Target (Partner 1)
    if not os.path.exists(CRBN_PDB):
        print(f"Error: Fixed partner PDB {CRBN_PDB} not found in current directory.")
        sys.exit(1)

    actpass_1 = "actpass_1"
    print(f"CRBN active residues: {crbn_active}")
    run_haddock3_restraints(CRBN_PDB, crbn_active, actpass_1)

    # 2. Generate actpass for CRBN (Partner 2 - Fixed)
    actpass_2 = "actpass_2"
    run_haddock3_restraints(args.target_pdb, args.target_res, actpass_2)

    # 3. Generate Ambiguous Restraints
    run_haddock3_ambig(actpass_1, actpass_2, args.output)

if __name__ == "__main__":
    main()
