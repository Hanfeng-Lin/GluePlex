import os
import glob
import subprocess
import sys
import yaml
import argparse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _gpu_pool import GPUPool, slot

# Template for the run script
SCRIPT_TEMPLATE = """import os
import subprocess
import sys
from pathlib import Path

YAML_CONTENT = \"\"\"
version: 1  # Optional, defaults to 1
sequences:
  - protein:
      id: A  # {protein_a_name}
      sequence: {protein_a_seq}
      msa: empty
  - protein:
      id: B  # {protein_b_name}
      sequence: {protein_b_seq}
      msa: empty
  - ligand:
      id: C
      smiles: '{ligand_smiles}'
templates:
    - cif: {cif_path} # if a pdb path is provided, Boltz will incrementally assign template chain ids based on the chain names in the PDB file (A1, A2, B1, etc)
      chain_id: [A, B]
      template_id: [A, B]
\"\"\"

def main(yaml_file, out_dir, output_format="pdb", diffusion_samples={diffusion_samples}, gpu_id=0, override=False):
    # the launcher may pin the device through the environment; only fall back
    # to the value baked in at generation time if it did not
    os.environ.setdefault("CUDA_VISIBLE_DEVICES", str(gpu_id))
    with open(yaml_file, "w") as f:
        f.write(YAML_CONTENT)
    print(f"Created {{yaml_file}}")

    cmd = [
        "boltz",
        "predict",
        yaml_file,
        "--out_dir",
        out_dir,
        "--output_format",
        output_format,
        "--diffusion_samples",
        str(diffusion_samples)
    ]

    # A fixed seed makes a run repeatable. Absent, boltz draws its own and the
    # pipeline behaves exactly as before, so adding this changes nothing for
    # configs that do not ask for it.
    if {boltz_seed} is not None:
        cmd += ["--seed", str({boltz_seed})]

    if override:
        cmd.append("--override")
    
    print(f"Running command: {{' '.join(cmd)}}")
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running boltz: {{e}}")
        sys.exit(1)

if __name__ == "__main__":
    main("{yaml_filename}", "{out_dir}", 
        output_format="pdb", diffusion_samples={diffusion_samples}, gpu_id={gpu_id}, override={override})
"""

def get_visible_devices(gpus=None):
    """Get list of visible CUDA devices."""
    if gpus:
        return gpus
        
    cuda_visible = os.environ.get("CUDA_VISIBLE_DEVICES")
    if cuda_visible:
        return [int(x.strip()) for x in cuda_visible.split(",")]
    else:
        # Fallback: try to count GPUs using nvidia-smi or just assume 0 if fails
        try:
            result = subprocess.run(["nvidia-smi", "-L"], capture_output=True, text=True)
            if result.returncode == 0:
                count = len(result.stdout.strip().split('\n'))
                return list(range(count))
        except FileNotFoundError:
            pass
        return [0] # Default to device 0 if detection fails

def process_cif(cif_file, gpu_id, cluster_num, config, output_dir, diffusion_samples, override=False, boltz_seed=None):
    """Generate script and run it."""
    cif_path = Path(cif_file)
    if not cif_path.exists():
        print(f"Warning: {cif_file} does not exist.")
        return

    # Use provided cluster_num
    cluster_name = f"cluster{cluster_num}"

    script_name = f"run_boltz2_{cluster_name}.py"
    
    protein_a_name = config.get('protein_a_name', 'Target')
    protein_b_name = config.get('protein_b_name', 'CRBN')
    
    yaml_filename = f"{protein_a_name}_{protein_b_name}_template{cluster_num}.yaml"
    
    # Output directory for this specific run
    out_dir = os.path.join(output_dir, f"{protein_a_name}_{protein_b_name}_template_noMSA_{diffusion_samples}runs")

    if override and os.path.exists(out_dir):
        import shutil
        print(f"Override set: Removing existing output directory {out_dir} to force re-processing.")
        shutil.rmtree(out_dir)

    script_content = SCRIPT_TEMPLATE.format(
        boltz_seed=boltz_seed,
        cif_path=f"{cif_path.resolve()}",
        yaml_filename=yaml_filename,
        out_dir=out_dir,
        gpu_id=gpu_id,
        protein_a_name=protein_a_name,
        protein_a_seq=config['protein_a_sequence'],
        protein_b_name=protein_b_name,
        protein_b_seq=config['protein_b_sequence'],
        ligand_smiles=config['ligand_smiles'],
        diffusion_samples=diffusion_samples,
        override=override
    )

    with open(script_name, "w") as f:
        f.write(script_content)
    
    print(f"Generated {script_name} for {cif_file} on GPU {gpu_id}")
    
    # Run the script
    cmd = ["python", script_name]
    try:
        subprocess.run(cmd, check=True)
        print(f"Finished {script_name}")
    except subprocess.CalledProcessError as e:
        print(f"Error running {script_name}: {e}")

def main():
    parser = argparse.ArgumentParser(description="Run Boltz batch processing.")
    parser.add_argument("--config", required=True, help="Path to project_config.yaml")
    parser.add_argument("--input_dir", required=True, help="Directory containing input CIF files.")
    parser.add_argument("--output_dir", required=True, help="Directory for output results.")
    parser.add_argument("--diffusion_samples", type=int, default=None, help="Number of diffusion samples.")
    parser.add_argument("--gpus", type=int, nargs="+", help="List of GPU IDs to use (e.g. 0 1).")
    parser.add_argument("--override", action="store_true", help="Override existing Boltz predictions.")
    
    args = parser.parse_args()
    
    # Load config
    with open(args.config, 'r') as f:
        full_config = yaml.safe_load(f)
        
    # Extract Step 3 config
    step3_config = full_config.get('step3', {})
    if not step3_config:
        print("Error: 'step3' section missing in config file.")
        sys.exit(1)

    # Use diffusion_samples from config if not overridden by CLI (CLI default is 20, so hard to tell if user passed it)
    # Let's assume if user passes --diffusion_samples it overrides config.
    # But argparse default makes it tricky.
    # Better: set argparse default to None.
    
    config_diffusion = step3_config.get('diffusion_samples', 20)
    if args.diffusion_samples is not None:
        diffusion_samples = args.diffusion_samples
    else:
        diffusion_samples = config_diffusion
        
    print(f"Using diffusion_samples: {diffusion_samples}")

    # Find CIF files
    cif_files = sorted(glob.glob(os.path.join(args.input_dir, "*.cif")))
    if not cif_files:
        print(f"No CIF files found in {args.input_dir}")
        return

    # Determine GPUs to use
    # Priority: CLI > Config > Auto-detect
    if args.gpus:
        gpus_list = args.gpus
    else:
        gpus_list = step3_config.get('gpus')
        
    gpus = get_visible_devices(gpus_list)
    print(f"Detected GPUs: {gpus}")
    
    # Hand each job whichever card is actually free rather than assigning by
    # i % len(gpus). Static round-robin both ignored other pipelines sharing the
    # machine and relied on mutating os.environ from several threads, which races
    # -- jobs piled onto one card while the rest idled.
    pool = GPUPool(gpus)

    def _run(cif_file, cluster_num):
        with slot(pool) as gpu_id:
            return process_cif(cif_file, gpu_id, cluster_num, step3_config,
                               args.output_dir, diffusion_samples, args.override,
                               boltz_seed=full_config.get("boltz_seed"))

    with ThreadPoolExecutor(max_workers=len(gpus)) as executor:
        futures = []
        for i, cif_file in enumerate(cif_files):
            # Use index + 1 as cluster number
            futures.append(executor.submit(_run, cif_file, i + 1))
        
        for future in futures:
            future.result()

if __name__ == "__main__":
    main()
