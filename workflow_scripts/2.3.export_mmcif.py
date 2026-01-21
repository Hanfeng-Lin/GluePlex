import os
import glob
import subprocess
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description="Export HADDOCK3 results to mmCIF.")
    parser.add_argument("--run_dir", default="run1-mdref-full", help="HADDOCK3 run directory.")
    
    args = parser.parse_args()
    
    # Define the directory and pattern
    # Assuming seletopclusts is always the folder of interest, but it might change index
    # Let's search for *seletopclusts* inside run_dir
    
    if not os.path.exists(args.run_dir):
        print(f"Error: Run directory {args.run_dir} does not exist.")
        sys.exit(1)
        
    # Find seletopclusts directory
    subdirs = glob.glob(os.path.join(args.run_dir, "*_seletopclusts"))
    if not subdirs:
        print(f"Error: No seletopclusts directory found in {args.run_dir}")
        sys.exit(1)
        
    # Use the first one found (usually there's only one)
    search_dir = subdirs[0]
    print(f"Using results directory: {search_dir}")
    
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
        sys.exit(1)
    except FileNotFoundError:
        print("Error: 'chimerax' executable not found in PATH.")
        sys.exit(1)

if __name__ == "__main__":
    main()
