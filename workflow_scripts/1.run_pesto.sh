#!/bin/bash
# Helper script to run PeSTo apply_model.py
#
# Usage: ./run_pesto.sh [OPTIONS]
#
# Options:
#   --pdb_id, -n ID       PDB ID to fetch (e.g., 1E9H or 1E9H_A)
#   --pdb_file, -f FILE   Path to input PDB file or directory
#   --output, -o DIR      Output directory (defaults to current directory)
#   --model_path, -m PATH Path to model directory
#   --device, -d DEVICE   Device to use (cpu or cuda)
#   --help, -h            Show help message
#
# Examples:
#   ./run_pesto.sh --pdb_id 1E9H_A
#   ./run_pesto.sh --pdb_file my_structure.pdb --output results/

# Get the directory of this script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PESTO_SCRIPT="$SCRIPT_DIR/../tools/PeSTo/apply_model.py"

# Run the python script with all arguments passed to this shell script
python3 "$PESTO_SCRIPT" "$@"
