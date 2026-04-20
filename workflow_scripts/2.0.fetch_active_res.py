import sys
import os
import subprocess

def run_naccess(pdb_file):
    """Runs naccess on the PDB file and returns the .rsa filename."""
    # naccess requires the file to be in the current directory or handled carefully
    # It outputs .rsa file in the current directory
    
    # We'll run naccess in the directory of the pdb_file to avoid path issues if possible,
    # or just run it where we are.
    # The user command was `naccess VAV_af2_pesto.pdb`
    
    abs_pdb = os.path.abspath(pdb_file)
    pdb_dir = os.path.dirname(abs_pdb)
    pdb_name = os.path.basename(abs_pdb)
    
    # Check if naccess is in path
    naccess_cmd = "naccess"
    
    # Run naccess
    # We run it in the pdb_dir so output files land there
    try:
        subprocess.run([naccess_cmd, pdb_name], cwd=pdb_dir, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except FileNotFoundError:
        print("Error: 'naccess' command not found. Please ensure it is in your PATH.")
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"Error running naccess: {e}")
        sys.exit(1)
        
    rsa_file = os.path.join(pdb_dir, os.path.splitext(pdb_name)[0] + ".rsa")
    return rsa_file

def get_sasa_residues(rsa_file):
    """
    Parses .rsa file and returns residues with > 15% accessibility.
    Logic:
    if (NF==13 && ($7>15 || $9>15))
    if (NF==14 && ($8>15 || $10>15))
    """
    sasa_res = set()
    if not os.path.exists(rsa_file):
        print(f"Error: RSA file {rsa_file} not found.")
        return sasa_res

    with open(rsa_file, 'r') as f:
        for line in f:
            if not line.startswith("RES"):
                continue
            
            parts = line.split()
            nf = len(parts)
            
            # parts[0] is RES
            # parts[1] is RES name (e.g. MET)
            # parts[2] is Chain or ResNum depending on format?
            # Let's rely on the user's column indices which are 1-based in awk, so 0-based in python:
            # $7 -> index 6, $9 -> index 8
            # $8 -> index 7, $10 -> index 9
            
            # We need to extract residue number.
            # Standard RSA format: RES ALA A 123 ...
            # NF=14 usually implies Chain is present: RES ALA A 123 ...
            # NF=13 usually implies Chain is missing: RES ALA 123 ...
            
            res_seq = None
            keep = False
            
            try:
                if nf == 13:
                    # RES ResName ResNum ...
                    # $3 is ResNum (index 2)
                    res_seq = parts[2]
                    val7 = float(parts[6]) # $7
                    val9 = float(parts[8]) # $9
                    if val7 > 15 or val9 > 15:
                        keep = True
                        
                elif nf == 14:
                    # RES ResName Chain ResNum ...
                    # $4 is ResNum (index 3)
                    res_seq = parts[3]
                    val8 = float(parts[7]) # $8
                    val10 = float(parts[9]) # $10
                    if val8 > 15 or val10 > 15:
                        keep = True
            except ValueError:
                continue
                
            if keep and res_seq:
                sasa_res.add(res_seq)
                
    return sasa_res

def get_pesto_residues(pdb_file):
    """
    Parses a PDB file and returns a set of residue indices (strings)
    for CA atoms with B-factor >= 0.3.
    """
    pesto_res = set()
    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith("ATOM"):
                    # Check for CA atom
                    atom_name = line[12:16].strip()
                    if atom_name == "CA":
                        try:
                            b_factor = float(line[60:66])
                            if b_factor >= 0.3:
                                res_seq = line[22:26].strip()
                                pesto_res.add(res_seq)
                        except ValueError:
                            continue
    except Exception as e:
        print(f"Error reading {pdb_file}: {e}")
        sys.exit(1)
        
    return pesto_res

def cleanup_naccess(rsa_file):
    """Removes naccess output files."""
    base = os.path.splitext(rsa_file)[0]
    for ext in ['.rsa', '.asa', '.log']:
        f = base + ext
        if os.path.exists(f):
            os.remove(f)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python 2.0.fetch_active_res.py <pdb_file>")
        sys.exit(1)
    
    pdb_file = sys.argv[1]
    
    # 1. Get PeSTo residues
    pesto_res = get_pesto_residues(pdb_file)
    # print(f"PeSTo residues: {len(pesto_res)}")
    
    # 2. Run NACCESS and get SASA residues
    sasa_res = set()
    try:
        rsa_file = run_naccess(pdb_file)
        sasa_res = get_sasa_residues(rsa_file)
    except SystemExit:
        print("Warning: naccess failed.", file=sys.stderr)

    # 3. Determine active residues with fallbacks
    if pesto_res and sasa_res:
        active_res = pesto_res.intersection(sasa_res)
        if not active_res:
            active_res = pesto_res
            print("Warning: PeSTo/SASA intersection empty, using PeSTo residues only.", file=sys.stderr)
    elif pesto_res:
        active_res = pesto_res
        print("Warning: No SASA data, using PeSTo residues only.", file=sys.stderr)
    elif sasa_res:
        active_res = sasa_res
        print("Warning: No PeSTo residues found, falling back to SASA surface residues.", file=sys.stderr)
    else:
        print("Error: Neither PeSTo nor SASA produced any residues.", file=sys.stderr)
        sys.exit(1)

    active_res = sorted(list(active_res), key=lambda x: int(x) if x.isdigit() else x)

    # Print as comma-separated list
    print(",".join(active_res))
