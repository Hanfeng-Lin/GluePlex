import os
import glob
import subprocess

def main():
    # Define the directory and pattern
    search_dir = "run1-mdref-full/10_seletopclusts"
    pattern = "cluster_*model_1.pdb.gz"
    search_path = os.path.join(search_dir, pattern)
    
    # Find all matching files
    files = glob.glob(search_path)
    
    if not files:
        print(f"No files found matching {search_path}")
        return

    print(f"Found {len(files)} files to convert.")
    
    # Create ChimeraX command file
    cxc_file = "convert_mmcif.cxc"
    with open(cxc_file, "w") as f:
        for pdb_file in files:
            # Construct output filename
            # e.g., run1-mdref-full/10_seletopclusts/cluster_1_model_1.pdb.gz -> ./cluster_1_model_1.cif
            base_name = os.path.basename(pdb_file).replace(".pdb.gz", "")
            output_file = f"{base_name}.pdb.cif"
            
            # Get absolute paths to be safe
            abs_pdb = os.path.abspath(pdb_file)
            abs_output = os.path.abspath(output_file)
            
            f.write(f"open {abs_pdb}\n")
            f.write(f"save {abs_output} format mmcif bestGuess true\n")
            f.write("close session\n")
        f.write("exit\n")
            
    print(f"Generated ChimeraX command file: {cxc_file}")
    
    # Run ChimeraX
    # Assuming chimerax is in the PATH
    command = ["chimerax", "--nogui", cxc_file]
    print(f"Running command: {' '.join(command)}")
    
    try:
        subprocess.run(command, check=True)
        print("Conversion completed successfully.")
    except subprocess.CalledProcessError as e:
        print("Error running ChimeraX:")
        print(e)
    except FileNotFoundError:
        print("Error: 'chimerax' executable not found in PATH.")

if __name__ == "__main__":
    main()
