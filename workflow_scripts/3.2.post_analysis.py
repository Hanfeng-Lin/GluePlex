import os
import re
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

def align_and_extract(results_dir, output_csv, align_chain="A", extract_chain="B"):
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
    
    headers = ["file", "cluster", "pass", "ipTM", "pair_chains_ipTM_avg"]
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

        m_pass2 = re.search(r"pass2_top(\d+)", basename)
        if m_pass2:
            cluster = f"p2_top{m_pass2.group(1)}"
            pass_num = 2
        else:
            m_t = re.search(r"template(\d+)", basename)
            cluster = m_t.group(1) if m_t else "0"
            pass_num = 1

        row = [basename, cluster, pass_num, iptm, pair_iptm] + extract_coords
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
    parser.add_argument("--align_chain", default="A", help="Chain ID to align on (Fixed Partner, e.g. CRBN).")
    parser.add_argument("--extract_chain", default="B", help="Chain ID to extract coordinates from (Mobile Target).")
    parser.add_argument("--ligand_name", default="Ligand", help="Name of the ligand.")
    parser.add_argument("--protein_name", default="Protein", help="Name of the target protein (protein B).")
    parser.add_argument("--protein_a_name", default="", help="Name of the fixed protein (protein A, e.g. CRBN).")
    parser.add_argument("--ref_pdb", default=None, help="Reference PDB for RMSD calculation.")
    parser.add_argument("--dim_reduction", default="umap", choices=["umap", "pca"],
                        help="Dimensionality reduction method: 'umap' (default) or 'pca'.")

    args = parser.parse_args()
    
    # 1. Collect results
    setup_results_dir(args.source_dir, args.results_dir)
    
    # 2. Align and extract
    align_and_extract(args.results_dir, args.output_csv, args.align_chain, args.extract_chain)

    # 2b. Compute reference RMSD (if --ref_pdb provided)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if args.ref_pdb:
        ref_script = os.path.join(script_dir, "_ref_align_rmsd.py")
        if os.path.exists(ref_script):
            print(f"Computing ref RMSD against {args.ref_pdb}...")
            cmd = [sys.executable, ref_script,
                   "--ref_pdb", args.ref_pdb,
                   "--input_csv", args.output_csv,
                   "--results_dir", args.results_dir,
                   "--align_chain", args.align_chain,
                   "--extract_chain", args.extract_chain]
            try:
                subprocess.run(cmd, check=True)
            except subprocess.CalledProcessError as e:
                print(f"Warning: ref RMSD computation failed: {e}")
        else:
            print(f"Warning: _ref_align_rmsd.py not found at {ref_script}")

    # 2c. Save top-10 log (after ref_rmsd so it can include that column)
    import pandas as _pd
    _df = _pd.read_csv(args.output_csv)
    _df_sorted = _df.sort_values("pair_chains_ipTM_avg", ascending=False).head(10)
    log_path = os.path.join(os.path.dirname(args.output_csv) or ".", "top10_models.log")
    with open(log_path, "w") as lf:
        has_rmsd = "ref_rmsd" in _df.columns
        hdr = f"{'Rank':<6}{'pair_ipTM':>12}{'ipTM':>10}"
        if has_rmsd:
            hdr += f"{'ref_rmsd':>10}"
        hdr += f"  {'File'}"
        lf.write(f"Top {len(_df_sorted)} models by pair_chains_ipTM_avg\n")
        lf.write(hdr + "\n")
        lf.write("-" * len(hdr) + "\n")
        for i, (_, r) in enumerate(_df_sorted.iterrows(), 1):
            line = f"{i:<6}{r['pair_chains_ipTM_avg']:>12.4f}{r['ipTM']:>10.4f}"
            if has_rmsd:
                line += f"{r['ref_rmsd']:>10.3f}" if not _pd.isna(r['ref_rmsd']) else f"{'NaN':>10}"
            line += f"  {r['file']}"
            lf.write(line + "\n")
    print(f"Top-10 model log saved to {log_path}")

    # 3. Run dimensionality reduction visualization
    method = args.dim_reduction.lower()
    script_map = {
        "umap": "_umap_coord_energy.py",
        "pca":  "_pca_coord_energy.py",
    }
    script_filename = script_map[method]
    print(f"Running {method.upper()} visualization ({script_filename})...")
    vis_script = os.path.join(script_dir, script_filename)

    if os.path.exists(vis_script):
        try:
            cmd = [sys.executable, vis_script,
                   "--ligand_name", args.ligand_name,
                   "--protein_name", args.protein_name,
                   "--protein_a_name", args.protein_a_name,
                   "--input_csv", args.output_csv]
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running {method.upper()} script: {e}")
    else:
        print(f"Warning: {method.upper()} script not found at {vis_script}")
    
    print("Post-analysis complete.")

if __name__ == "__main__":
    main()
