import os
import subprocess
import shutil
import glob
import json
import csv
import gemmi
import argparse
import sys

def setup_results_dir(source_dir, results_dir):
    """Collect PDB and JSON files into results directory."""
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    
    print(f"Collecting files from {source_dir} to {results_dir}...")
    
    # Copy PDBs recursively
    # Pattern: source_dir/**/predictions/**/*.pdb
    pdb_files = glob.glob(os.path.join(source_dir, "**", "*.pdb"), recursive=True)
    if not pdb_files:
        print(f"Warning: No PDB files found in {source_dir}")
    
    for pdb_file in pdb_files:
        shutil.copy(pdb_file, results_dir)
        
    # Copy JSONs recursively
    json_files = glob.glob(os.path.join(source_dir, "**", "*.json"), recursive=True)
    if not json_files:
        print(f"Warning: No JSON files found in {source_dir}")

    for json_file in json_files:
        shutil.copy(json_file, results_dir)

def get_ca_coords(structure, chain_id):
    """Extract CA coordinates for a specific chain."""
    coords = []
    model = structure[0]
    for chain in model:
        if chain.name == chain_id:
            for res in chain:
                ca = res.find_atom("CA", '*')
                if ca:
                    coords.extend([ca.pos.x, ca.pos.y, ca.pos.z])
            break
    return coords

def align_and_extract(results_dir, output_csv, align_chain="B", extract_chain="A"):
    """Align structures and extract coordinates/scores."""
    pdb_files = sorted(glob.glob(os.path.join(results_dir, "*.pdb")))
    if not pdb_files:
        print("No PDB files found.")
        return

    ref_pdb = pdb_files[0]
    print(f"Using {os.path.basename(ref_pdb)} as reference for alignment.")
    
    ref_st = gemmi.read_structure(ref_pdb)
    
    # Prepare CSV header
    # We need to know the number of CA atoms in extract_chain to generate headers
    ref_extract_coords = get_ca_coords(ref_st, extract_chain)
    print(f"Reference Chain {extract_chain} coords length: {len(ref_extract_coords)}")
    
    num_atoms = len(ref_extract_coords) // 3
    
    headers = ["file", "cluster", "ipTM", "pair_chains_ipTM_CRBN_VAV1_avg"]
    for i in range(1, num_atoms + 1):
        headers.extend([f"{i}x", f"{i}y", f"{i}z"])
        
    data_rows = []
    
    for pdb_file in pdb_files:
        basename = os.path.basename(pdb_file)
        # print(f"Processing {basename}...")
        
        # Load structure
        st = gemmi.read_structure(pdb_file)
        
        # Align to reference align_chain
        try:
            # Perform superposition using Chain objects
            # Use CaP selection (CA for proteins, P for nucleic acids)
            sup = gemmi.calculate_superposition(ref_st[0][align_chain].whole(), st[0][align_chain].whole(), gemmi.PolymerType.PeptideL, gemmi.SupSelect.CaP)
            st[0].transform_pos_and_adp(sup.transform)
        except Exception as e:
            print(f"Alignment failed for {basename}: {e}")
            continue
            
        # Extract extract_chain CA coords
        extract_coords = get_ca_coords(st, extract_chain)
        
        if len(extract_coords) != len(ref_extract_coords):
            print(f"Warning: Atom count mismatch for {basename}. Ref {extract_chain}: {len(ref_extract_coords)}, Curr {extract_chain}: {len(extract_coords)}. Skipping.")
            continue
            
        # Load scores
        # The JSON file has the same basename but with .json extension and 'confidence_' prefix
        # e.g. CRBN_vav1_template1_model_0.pdb -> confidence_CRBN_vav1_template1_model_0.json
        # Adjusting prefix logic to be more generic if possible, or assume standard Boltz output
        # Boltz output: confidence_{name}.json
        
        # Try to find matching json
        # name is usually {proteinA}_{proteinB}_template{N}_model_{M}
        name_no_ext = os.path.splitext(basename)[0]
        json_basename = f"confidence_{name_no_ext}.json"
        json_file = os.path.join(results_dir, json_basename)
        
        iptm = 0.0
        pair_iptm = 0.0
        
        if os.path.exists(json_file):
            with open(json_file, 'r') as f:
                scores = json.load(f)
                iptm = scores.get("iptm", 0.0)
                # Calculate average of pair_chains_iptm 0->1 and 1->0
                pair_chains = scores.get("pair_chains_iptm", {})
                val_0_1 = pair_chains.get("0", {}).get("1", 0.0)
                val_1_0 = pair_chains.get("1", {}).get("0", 0.0)
                pair_iptm = (val_0_1 + val_1_0) / 2.0
        else:
            # Fallback for older naming convention if needed, but standard boltz is consistent
            print(f"Warning: JSON score file not found for {basename} at {json_file}")

        try:
            cluster = basename.split("template")[1].split("_")[0]
        except:
            cluster = "0"

        row = [basename, cluster, iptm, pair_iptm] + extract_coords
        data_rows.append(row)
        
    # Write CSV
    print(f"Writing {len(data_rows)} rows to {output_csv}")
    with open(output_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        writer.writerows(data_rows)
        
    print(f"Saved aligned coordinates and scores to {output_csv}")

def main():
    parser = argparse.ArgumentParser(description="Run Boltz post-analysis.")
    parser.add_argument("--source_dir", required=True, help="Directory containing Boltz output folders.")
    parser.add_argument("--results_dir", default="results", help="Directory to collect flattened results.")
    parser.add_argument("--output_csv", default="atom.csv", help="Output CSV filename.")
    parser.add_argument("--align_chain", default="B", help="Chain ID to align on (Fixed Partner).")
    parser.add_argument("--extract_chain", default="A", help="Chain ID to extract coordinates from (Mobile Target).")
    parser.add_argument("--ligand_name", default="Ligand", help="Name of the ligand.")
    parser.add_argument("--protein_name", default="Protein", help="Name of the protein.")
    
    args = parser.parse_args()
    
    # 1. Collect results
    setup_results_dir(args.source_dir, args.results_dir)
    
    # 2. Align and extract
    align_and_extract(args.results_dir, args.output_csv, args.align_chain, args.extract_chain)
    
    # 3. Run UMAP visualization
    print("Running UMAP visualization...")
    # Assuming _umap_coord_energy.py is in the same directory as this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    umap_script = os.path.join(script_dir, "_umap_coord_energy.py")
    
    if os.path.exists(umap_script):
        try:
            cmd = [sys.executable, umap_script, "--ligand_name", args.ligand_name, "--protein_name", args.protein_name]
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running UMAP script: {e}")
    else:
        print(f"Warning: UMAP script not found at {umap_script}")
    
    print("Post-analysis complete.")

if __name__ == "__main__":
    main()
