import subprocess
import os
import argparse
import sys

# ==========================================
# Configuration Parameters
# ==========================================
# Fixed parameters for CRBN (Partner 1)
CRBN_PDB = "CRBN_pesto.pdb"
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
    parser.add_argument("--target_res", default="", help="Comma-separated active residues for the target.")
    parser.add_argument("--target_unrestrained", action="store_true", help="Leave target unrestrained; only CRBN side has restraints.")
    parser.add_argument("--output", default="restraint.tbl", help="Output filename for the restraints table.")
    
    args = parser.parse_args()
    
    # 1. Generate actpass for CRBN (always restrained)
    if not os.path.exists(CRBN_PDB):
        print(f"Error: Fixed partner PDB {CRBN_PDB} not found in current directory.")
        sys.exit(1)

    actpass_1 = "actpass_1"
    run_haddock3_restraints(CRBN_PDB, CRBN_ACTIVE_RES, actpass_1)

    # 2. Generate actpass for Target
    actpass_2 = "actpass_2"
    if args.target_unrestrained:
        # Write an empty actpass file (no active, no passive) so CRBN drives the docking
        print("Target is unrestrained — writing empty actpass for target side.")
        with open(actpass_2, "w") as f:
            f.write("\n\n")
    else:
        run_haddock3_restraints(args.target_pdb, args.target_res, actpass_2)

    # 3. Generate Ambiguous Restraints
    run_haddock3_ambig(actpass_1, actpass_2, args.output)

if __name__ == "__main__":
    main()
