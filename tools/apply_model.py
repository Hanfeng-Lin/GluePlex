#!/usr/bin/env python3
import os
import sys
import argparse
import urllib.request
import torch as pt
from tqdm import tqdm
import re
from glob import glob

# Determine the directory of the script
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Add the script directory to sys.path to allow importing src
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

# Add src to path if needed, though usually running from root works if src is there.
# The notebook imports from src.dataset etc, so we assume we are running from the repo root.

from src.dataset import StructuresDataset, collate_batch_features
from src.data_encoding import encode_structure, encode_features, extract_topology
from src.structure import encode_bfactor, concatenate_chains, split_by_chain
from src.structure_io import save_pdb, read_pdb

def load_model(model_path, device):
    # Add model path to sys.path to allow importing config, data_handler, model
    if model_path not in sys.path:
        sys.path.insert(0, model_path)
    
    try:
        from config import config_model
        from model import Model
    except ImportError as e:
        print(f"Error importing model modules from {model_path}: {e}")
        sys.exit(1)

    # create model
    model = Model(config_model)

    # reload model
    model_ckpt = os.path.join(model_path, 'model_ckpt.pt')
    if not os.path.exists(model_ckpt):
         # Fallback to model.pt if model_ckpt.pt doesn't exist
         model_ckpt = os.path.join(model_path, 'model.pt')
    
    if not os.path.exists(model_ckpt):
        print(f"Model checkpoint not found at {model_ckpt}")
        sys.exit(1)

    print(f"Loading model from {model_ckpt}...")
    model.load_state_dict(pt.load(model_ckpt, map_location=pt.device("cpu")))

    # set model to inference
    model = model.eval().to(device)
    return model

def fetch_pdb(pdb_id, output_dir):
    # parse pdb id
    if "_" in pdb_id:
        pdb_code, chain_id = pdb_id.split("_")
    else:
        pdb_code = pdb_id
        chain_id = None

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # paths
    pdb_filepath = os.path.join(output_dir, "{}.pdb".format(pdb_code))
    final_filepath = os.path.join(output_dir, "{}.pdb".format(pdb_id))

    # download
    url = "https://files.rcsb.org/download/{}.pdb".format(pdb_code)
    try:
        print(f"Downloading {pdb_code} from {url}...")
        urllib.request.urlretrieve(url, pdb_filepath)
        print(f"Downloaded {pdb_code} to {pdb_filepath}")
    except Exception as e:
        print(f"Failed to download PDB {pdb_code}: {e}")
        return None

    # extract chain if needed
    if chain_id is not None:
        # read structure
        structure = read_pdb(pdb_filepath)
        
        # split chains
        chains = split_by_chain(structure)
        
        # search for chain (handling model suffix like A:0)
        target_chain_key = None
        for key in chains.keys():
            # key format is usually CHAIN:MODEL
            if key.split(":")[0] == chain_id:
                target_chain_key = key
                break
        
        # check if chain exists
        if target_chain_key:
            # save specific chain
            save_pdb({target_chain_key: chains[target_chain_key]}, final_filepath)
            print(f"Extracted chain {chain_id} (found as {target_chain_key}) to {final_filepath}")
            
            # remove original file
            os.remove(pdb_filepath)
            return final_filepath
        else:
            print(f"Chain {chain_id} not found in {pdb_code}. Available chains: {list(chains.keys())}")
            return None
    else:
        # If no chain specified, the downloaded file is the final file
        # But if pdb_id was just code, pdb_filepath is final_filepath
        return pdb_filepath

def process_data(data_path, model, device, output_path=None):
    # If data_path is a list (from fetching), use it directly
    if isinstance(data_path, list):
        pdb_filepaths = data_path
    else:
        # find pdb files and ignore already predicted ones
        if os.path.isdir(data_path):
            pdb_filepaths = glob(os.path.join(data_path, "*.pdb"), recursive=True)
        else:
            pdb_filepaths = [data_path]

    # Skip PeSTo's own per-interface outputs (..._i1.pdb, _i2.pdb, ...). Match the
    # BASENAME with an explicit pattern: testing `"_i" not in fp` against the whole
    # path silently discards every input whose *directory* happens to contain "_i"
    # (e.g. an output folder named "runs_improve"), leaving "No PDB files found".
    pdb_filepaths = [fp for fp in pdb_filepaths
                     if not re.search(r"_i\d+\.pdb$", os.path.basename(fp))]
    
    if not pdb_filepaths:
        print(f"No PDB files found to process.")
        return

    print(f"Found {len(pdb_filepaths)} PDB files to process.")

    # create dataset loader with preprocessing
    dataset = StructuresDataset(pdb_filepaths, with_preprocessing=True)

    # run model on all subunits
    with pt.no_grad():
        for subunits, filepath in tqdm(dataset):
            # concatenate all chains together
            structure = concatenate_chains(subunits)

            # encode structure and features
            X, M = encode_structure(structure)
            #q = pt.cat(encode_features(structure), dim=1)
            q = encode_features(structure)[0]

            # extract topology
            ids_topk, _, _, _, _ = extract_topology(X, 64)

            # pack data and setup sink (IMPORTANT)
            X, ids_topk, q, M = collate_batch_features([[X, ids_topk, q, M]])

            # run model
            z = model(X.to(device), ids_topk.to(device), q.to(device), M.float().to(device))

            # for all predictions
            for i in range(z.shape[1]):
                # prediction
                p = pt.sigmoid(z[:,i])

                # encode result
                structure = encode_bfactor(structure, p.cpu().numpy())

                # save results
                # Construct output filename
                base_name = os.path.basename(filepath)
                name_no_ext = os.path.splitext(base_name)[0]
                
                if i == 0:
                    output_filename = f"{name_no_ext}_pesto.pdb"
                else:
                    output_filename = f"{name_no_ext}_i{i}.pdb"
                
                # Determine final output directory
                if output_path:
                    final_out_dir = output_path
                else:
                    final_out_dir = os.path.dirname(filepath)
                
                if not os.path.exists(final_out_dir):
                    os.makedirs(final_out_dir)

                final_output_path = os.path.join(final_out_dir, output_filename)
                save_pdb(split_by_chain(structure), final_output_path)

def main():
    parser = argparse.ArgumentParser(description="Apply PeSTo model to PDB files.")
    parser.add_argument("--pdb_id", "-n", type=str, help="PDB ID to fetch (e.g., 1E9H or 1E9H_A).")
    parser.add_argument("--pdb_file", "-f", type=str, help="Path to input PDB file or directory.")
    parser.add_argument("--output", "-o", type=str, default=".", help="Path to output directory. Defaults to current directory.")
    parser.add_argument("--model_path", "-m", type=str, default="model/save/i_v4_1_2021-09-07_11-21", help="Path to the model directory.")
    parser.add_argument("--device", "-d", type=str, default="cpu", help="Device to use (cpu or cuda).")

    args = parser.parse_args()

    if not args.pdb_id and not args.pdb_file:
        parser.error("Either --pdb_id or --pdb_file must be provided.")

    device = pt.device(args.device)
    
    # Check if model path exists
    # If model path is relative, check relative to SCRIPT_DIR first
    if not os.path.isabs(args.model_path):
        abs_model_path = os.path.join(SCRIPT_DIR, args.model_path)
        if os.path.exists(abs_model_path):
            args.model_path = abs_model_path
    
    if not os.path.exists(args.model_path):
        print(f"Model path {args.model_path} does not exist. Please check the path.")
        sys.exit(1)

    model = load_model(args.model_path, device)

    pdb_files_to_process = []
    
    # Determine output directory
    # If output is not specified, we need a default.
    # For pdb_id, default could be current dir.
    # For pdb_file, default is dirname of file.
    # process_data handles output_path=None by using dirname of input file.
    # But for pdb_id, we download somewhere.
    
    output_dir = args.output
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if args.pdb_id:
        # Download to output_dir if specified, else current directory?
        # Let's say if output_dir is None, download to current dir.
        download_dir = output_dir if output_dir else "."
        fetched_pdb = fetch_pdb(args.pdb_id, download_dir)
        if fetched_pdb:
            pdb_files_to_process.append(fetched_pdb)
    
    if args.pdb_file:
        if os.path.isdir(args.pdb_file):
             pdb_files_to_process.extend(glob(os.path.join(args.pdb_file, "*.pdb"), recursive=True))
        elif os.path.isfile(args.pdb_file):
             pdb_files_to_process.append(args.pdb_file)
        else:
            print(f"File or directory {args.pdb_file} not found.")

    if pdb_files_to_process:
        # If we fetched a PDB, we passed it as a list.
        # process_data expects a path or list.
        process_data(pdb_files_to_process, model, device, output_dir)

if __name__ == "__main__":
    main()
