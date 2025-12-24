import subprocess
import os

# ==========================================
# Configuration Parameters
# Use E3 as chain A
# ==========================================
TASKS = [
    {
        "pdb_file": "CRBN_pesto.pdb",
        "active_residues": "105,332,341",
        "output_filename": "actpass_1"
    },
    {
        "pdb_file": "6S76_pesto.pdb",
        "active_residues": "3,4,5,6,7,8,9,35,36,37,38,60,64,67,68,69,70,71,72,76,77,78,79,80,81,82,83,92,134,135,136",
        "output_filename": "actpass_2"
    },
]
# ==========================================

def run_haddock3_restraints(pdb_file, active_residues, output_filename):
    """
    Runs haddock3-restraints passive_from_active and writes the output to a file.
    
    Args:
        pdb_file (str): Path to the PDB file.
        active_residues (str): Comma-separated string of active residues.
        output_filename (str): Name of the output file to write.
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
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

def run_haddock3_ambig(actpass_file1, actpass_file2, output_filename):
    """
    Runs haddock3-restraints active_passive_to_ambig and writes the output to a file.
    
    Args:
        actpass_file1 (str): Path to the first actpass file.
        actpass_file2 (str): Path to the second actpass file.
        output_filename (str): Name of the output file to write.
    """
    command = ["/home/hlin/miniconda3/envs/haddock3/bin/haddock3-restraints", "active_passive_to_ambig", actpass_file1, actpass_file2]
    
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
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

def main():
    generated_files = []
    for task in TASKS:
        run_haddock3_restraints(
            pdb_file=task["pdb_file"],
            active_residues=task["active_residues"],
            output_filename=task["output_filename"]
        )
        generated_files.append(task["output_filename"])

    if len(generated_files) == 2:
        run_haddock3_ambig(generated_files[0], generated_files[1], "restraint.tbl")

if __name__ == "__main__":
    main()
