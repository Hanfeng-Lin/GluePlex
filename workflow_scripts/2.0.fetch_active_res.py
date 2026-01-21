import sys
import os

def get_active_residues(pdb_file):
    """
    Parses a PDB file and returns a list of residue indices (strings)
    for CA atoms with B-factor >= 0.3.
    """
    active_res = []
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
                                active_res.append(res_seq)
                        except ValueError:
                            continue
    except FileNotFoundError:
        print(f"Error: File {pdb_file} not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading {pdb_file}: {e}")
        sys.exit(1)
        
    return active_res

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python 0.fetch_active_res.py <pdb_file>")
        sys.exit(1)
    
    pdb_file = sys.argv[1]
    active_residues = get_active_residues(pdb_file)
    
    # Print as comma-separated list
    print(",".join(active_residues))
