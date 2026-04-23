import os
import sys
import yaml
import argparse
import subprocess
from pathlib import Path

def load_config(config_path):
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def run_command(cmd, env=None, cwd=None, log_file=None):
    """Runs a command in a specific conda environment and logs output."""
    command = []
    if env:
        command = ["conda", "run", "-n", env, "--no-capture-output"]
    command.extend(cmd)
    
    cmd_str = " ".join(command)
    print(f"Running: {cmd_str}")
    
    if log_file:
        print(f"Logging to: {log_file}")
        with open(log_file, "a") as f:
            f.write(f"\n{'='*40}\n")
            f.write(f"Command: {cmd_str}\n")
            f.write(f"CWD: {cwd if cwd else os.getcwd()}\n")
            f.write(f"{'='*40}\n")
            f.flush()
            
            try:
                # Run and capture output, writing to both stdout and log file
                process = subprocess.Popen(
                    command, 
                    cwd=cwd, 
                    stdout=subprocess.PIPE, 
                    stderr=subprocess.STDOUT, # Merge stderr into stdout
                    text=True,
                    bufsize=1
                )
                
                for line in process.stdout:
                    print(line, end='') # Print to console
                    f.write(line)       # Write to log
                
                process.wait()
                if process.returncode != 0:
                    f.write(f"\nCommand failed with return code {process.returncode}\n")
                    print(f"Error: Command failed with return code {process.returncode}")
                    sys.exit(1)
                    
            except Exception as e:
                f.write(f"\nException: {e}\n")
                print(f"Error running command: {e}")
                sys.exit(1)
    else:
        # Legacy behavior (no log file)
        try:
            subprocess.run(command, check=True, cwd=cwd)
        except subprocess.CalledProcessError as e:
            print(f"Error running command: {e}")
            sys.exit(1)

def step1_pesto(config, project_dir, input_dir, log_file):
    # Clear log file
    if log_file:
        with open(log_file, 'w') as f:
            f.write("")
            
    print("\n=== Step 1: PeSTo ===")
    target_pdb = config['target_pdb']
    
    # Output directory for Step 1
    step1_dir = os.path.join(project_dir, "step1_pesto")
    if not os.path.exists(step1_dir):
        os.makedirs(step1_dir)
    
    # Determine if it's an ID or file
    if os.path.exists(target_pdb):
        # It's a file. Copy it to input_dir.
        import shutil
        dest_pdb = os.path.join(input_dir, os.path.basename(target_pdb))
        if os.path.abspath(target_pdb) != os.path.abspath(dest_pdb):
            shutil.copy(target_pdb, dest_pdb)
            print(f"Copied input PDB to {dest_pdb}")
        
        # Use the copied file for processing
        arg = ["--pdb_file", dest_pdb]
        base = os.path.splitext(os.path.basename(dest_pdb))[0]
    else:
        # It's an ID. Fetch it.
        # We run apply_model.py with cwd=input_dir so it fetches there
        arg = ["--pdb_id", target_pdb]
        base = target_pdb
        
    # Output to step1_dir
    # We need to use absolute path for apply_model.py since we might change cwd
    script_path = os.path.abspath("tools/PeSTo/apply_model.py")
    
    cmd = ["python", script_path] + arg + ["--output", step1_dir]
    
    # Run in input_dir so fetched PDBs land there
    run_command(cmd, env="pesto", cwd=input_dir, log_file=log_file)
    
    # Determine output filename
    # If ID: {ID}_pesto.pdb
    # If File: {basename}_pesto.pdb
    output_pdb = os.path.join(step1_dir, f"{base}_pesto.pdb")
    if not os.path.exists(output_pdb):
        print(f"Error: Expected output {output_pdb} not found.")
        sys.exit(1)
        
    print(f"Step 1 Complete. Output: {output_pdb}")
    return output_pdb

def step2_haddock(config, target_pesto_pdb, project_dir, log_file):
    # Clear log file
    if log_file:
        with open(log_file, 'w') as f:
            f.write("")
            
    print("\n=== Step 2: HADDOCK3 ===")
    
    # Output directory for Step 2
    step2_dir = os.path.join(project_dir, "step2_haddock")
    if not os.path.exists(step2_dir):
        os.makedirs(step2_dir)
    
    # 2.0 Fetch active residues (Optional/Placeholder logic)
    active_res = config.get('target_active_res')
    if not active_res:
        print("Fetching active residues using 2.0.fetch_active_res.py...")
        # Run in pesto env
        cmd = ["conda", "run", "-n", "pesto", "python", "workflow_scripts/2.0.fetch_active_res.py", target_pesto_pdb]
        
        # We need to capture output for this specific command to get the residues
        # run_command doesn't return output, so we use subprocess directly here
        # But we should log it too
        cmd_str = " ".join(cmd)
        print(f"Running: {cmd_str}")
        if log_file:
            with open(log_file, "a") as f:
                f.write(f"\n{'='*40}\nCommand: {cmd_str}\n{'='*40}\n")
        
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        active_res = result.stdout.strip()
        
        if log_file:
            with open(log_file, "a") as f:
                f.write(result.stdout)
                f.write(result.stderr)
                
        print(f"Fetched active residues: {active_res}")
    
    if not active_res:
        print("Error: No active residues found or provided.")
        sys.exit(1)

    # Copy CRBN_pesto.pdb to step2_dir
    # Assuming it's in the base directory (current working directory of the script)
    crbn_pdb = "CRBN_pesto.pdb"
    if os.path.exists(crbn_pdb):
        import shutil
        shutil.copy(crbn_pdb, step2_dir)
        print(f"Copied {crbn_pdb} to {step2_dir}")
    else:
        print(f"Warning: {crbn_pdb} not found in base directory. Step 2 might fail if it's not in {step2_dir}.")

    # 2.1 Generate Restraints
    # We run this inside step2_dir to keep generated files there
    print("Generating restraints...")
    restraint_file = "restraint.tbl"

    # Need absolute path for target_pesto_pdb if changing cwd
    abs_target_pesto_pdb = os.path.abspath(target_pesto_pdb)

    # Script path needs to be absolute or relative to step2_dir
    # Assuming run_glueplex.py is in base_dir
    script_path = os.path.abspath("workflow_scripts/2.1.generate_actpass.py")

    cmd = ["python", script_path,
           "--target_pdb", abs_target_pesto_pdb,
           "--target_res", active_res]
    cmd.extend(["--output", restraint_file])

    run_command(cmd, env="haddock3", cwd=step2_dir, log_file=log_file)
    
    # 2.2 Run HADDOCK3 Workflow
    print("Running HADDOCK3 docking...")
    
    # HADDOCK run directory inside step2_dir
    run_dir_name = config.get('haddock_run_dir', 'haddock_results')
    # If run_dir_name is just a name, put it in step2_dir
    # If it's a path, respect it? User wants everything in parent folder.
    # Let's force it into step2_dir
    run_dir = os.path.join(step2_dir, os.path.basename(run_dir_name))
    
    ncores = config.get('haddock_ncores', 48)
    script_path_22 = os.path.abspath("workflow_scripts/2.2.run_haddock_workflow.py")
    
    cmd = ["python", script_path_22,
           "--target_pdb", abs_target_pesto_pdb,
           "--restraint", restraint_file,
           "--run_dir", run_dir,
           "--ncores", str(ncores)]
           
    run_command(cmd, env="haddock3", cwd=step2_dir, log_file=log_file)
    
    # 2.3 Export mmCIF
    print("Exporting results to mmCIF...")
    script_path_23 = os.path.abspath("workflow_scripts/2.3.export_mmcif.py")
    cmd = ["python", script_path_23, "--run_dir", run_dir]
    run_command(cmd, env="haddock3", cwd=step2_dir, log_file=log_file) 
    
    print("Step 2 Complete.")
    return run_dir

def step3_boltz(config, haddock_run_dir, config_path, project_dir, log_file, override=False, ref_pdb=None):
    # Clear log file
    if log_file:
        with open(log_file, 'w') as f:
            f.write("")
            
    print("\n=== Step 3: Boltz ===")
    
    # Output directory for Step 3
    step3_dir = os.path.join(project_dir, "step3_boltz")
    if not os.path.exists(step3_dir):
        os.makedirs(step3_dir)
    
    # 3.1 Run Batch Clusters
    # Input dir for CIF files is step2_haddock (where 2.3 exported them)
    step2_dir = os.path.join(project_dir, "step2_haddock")
    if not os.path.exists(step2_dir):
         print(f"Error: Step 2 directory {step2_dir} not found.")
         sys.exit(1)
    input_dir = step2_dir
    
    # Boltz output dir inside step3_dir
    # config['boltz_output_dir'] might be "boltz_results"
    boltz_out_name = config.get('boltz_output_dir', 'boltz_results')
    output_dir = os.path.join(step3_dir, os.path.basename(boltz_out_name))
    
    print(f"Running Boltz batch processing from {input_dir}...")
    
    script_path_31 = os.path.abspath("workflow_scripts/3.1.run_batch_clusters.py")
    abs_config_path = os.path.abspath(config_path)
    
    cmd = ["python", script_path_31,
           "--config", abs_config_path,
           "--input_dir", input_dir,
           "--output_dir", output_dir]
           
    if override:
        cmd.append("--override")
    
    run_command(cmd, env="boltz_analysis", cwd=step3_dir, log_file=log_file)
    
    # 3.2 Post Analysis
    print("Running Post Analysis...")
    
    # Construct final results directory dynamically
    step3_config = config.get('step3', {})
    protein_a = step3_config.get('protein_a_name', 'Target')
    protein_b = step3_config.get('protein_b_name', 'CRBN')
    
    # Read ligand name
    ligand_name = step3_config.get('ligand_name', 'Ligand')
    
    # Final results inside step3_dir
    final_results_dir = os.path.join(step3_dir, f"{protein_a}-{protein_b}")
    
    print(f"Final results will be saved to: {final_results_dir}")
    
    script_path_32 = os.path.abspath("workflow_scripts/3.2.post_analysis.py")
    
    cmd = ["python", script_path_32,
           "--source_dir", output_dir,
           "--results_dir", final_results_dir,
           "--align_chain", "A",
           "--extract_chain", "B",
           "--ligand_name", ligand_name,
           "--protein_name", protein_b,
           "--protein_a_name", protein_a,
           "--dim_reduction", "all"]
    if ref_pdb:
        cmd.extend(["--ref_pdb", ref_pdb])
    run_command(cmd, env="boltz_analysis", cwd=step3_dir, log_file=log_file)

    print("Step 3 Complete.")

def step4_double_pass(config, config_path, project_dir, log_file, override=False, ref_pdb=None):
    """Pass-2 Boltz refinement: top-N from step-3 atom.csv as templates, then merged post-analysis."""
    if log_file:
        with open(log_file, 'w') as f:
            f.write("")

    print("\n=== Step 4: Double-Pass Boltz ===")

    step3_dir = os.path.join(project_dir, "step3_boltz")
    step4_dir = os.path.join(project_dir, "step4_boltz")
    if not os.path.exists(step3_dir):
        print(f"Error: step 3 directory {step3_dir} not found — run step 3 first.")
        sys.exit(1)
    if not os.path.exists(step4_dir):
        os.makedirs(step4_dir)

    abs_config_path = os.path.abspath(config_path)

    # 4.1 Prepare pass-2 jobs
    print("Preparing pass-2 jobs (4.1)...")
    script_41 = os.path.abspath("workflow_scripts/4.1.prepare_second_pass.py")
    cmd = ["python", script_41,
           "--config", abs_config_path,
           "--step3_dir", step3_dir,
           "--step4_dir", step4_dir]
    if override:
        cmd.append("--override")
    run_command(cmd, env="boltz_analysis", cwd=step4_dir, log_file=log_file)

    # 4.2 Run pass-2 Boltz jobs
    print("Running pass-2 Boltz jobs (4.2)...")
    script_42 = os.path.abspath("workflow_scripts/4.2.run_second_pass.py")
    cmd = ["python", script_42, "--step4_dir", step4_dir]
    run_command(cmd, env="boltz_analysis", cwd=step4_dir, log_file=log_file)

    # 4.3 Merge + unified post-analysis
    print("Merging pass-1 + pass-2 and running post-analysis (4.3)...")
    script_43 = os.path.abspath("workflow_scripts/4.3.merge_and_analyze.py")
    cmd = ["python", script_43,
           "--config", abs_config_path,
           "--step3_dir", step3_dir,
           "--step4_dir", step4_dir,
           "--output_csv", "atom_2ndpass.csv",
           "--dim_reduction", "all"]
    if ref_pdb:
        cmd.extend(["--ref_pdb", ref_pdb])
    run_command(cmd, env="boltz_analysis", cwd=step4_dir, log_file=log_file)

    print("Step 4 Complete.")


def main():
    parser = argparse.ArgumentParser(description="GluePlex Workflow Orchestrator")
    parser.add_argument("--config", default="project_config.yaml", help="Path to project configuration file.")
    parser.add_argument("--step", type=int, default=0, help="Step to run (1, 2, 3, 4). 0 runs all.")
    parser.add_argument("--override", action="store_true", help="Override existing Boltz predictions.")
    parser.add_argument("--doublepass", default=True, action=argparse.BooleanOptionalAction,
                        help="Run step 4 (double-pass Boltz refinement) after step 3. "
                             "Enabled by default; use --no-doublepass to skip.")
    parser.add_argument("--ref", default=None,
                        help="Reference PDB for RMSD calculation. Must have 3 chains (A, B, ligand).")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.config):
        print(f"Error: Config file {args.config} not found.")
        sys.exit(1)
        
    config = load_config(args.config)

    # Validate reference PDB if provided
    abs_ref = None
    if args.ref:
        if not os.path.exists(args.ref):
            print(f"Error: Reference PDB {args.ref} not found.")
            sys.exit(1)
        abs_ref = os.path.abspath(args.ref)
        chains = set()
        with open(abs_ref) as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")) and len(line) > 21:
                    chains.add(line[21])
        if len(chains) != 3:
            print(f"Error: Reference PDB must have exactly 3 chains, found {len(chains)}: {sorted(chains)}")
            sys.exit(1)
        print(f"Reference PDB: {abs_ref} (chains: {sorted(chains)})")

    # Create Project Directory.
    # On the doublepass branch, the canonical layout is
    # Result_{protein_a_name}_{protein_b_name}_GluePlex2pass — derive it from
    # step3 names so configs don't need to hardcode the project name.
    step3_cfg = config.get('step3', {}) or {}
    a_name = step3_cfg.get('protein_a_name')
    b_name = step3_cfg.get('protein_b_name')
    lig_name = step3_cfg.get('ligand_name')
    if a_name and b_name:
        # Convention on disk: target (b) first, then E3 (a), then ligand
        if lig_name:
            project_name = f"Result_{b_name}_{a_name}_{lig_name}_GluePlex2pass"
        else:
            project_name = f"Result_{b_name}_{a_name}_GluePlex2pass"
    else:
        project_name = config.get('project_name', 'GluePlex_Project')
    base_dir = config.get('base_dir', os.getcwd())
    project_dir = os.path.join(base_dir, project_name)
    
    if not os.path.exists(project_dir):
        os.makedirs(project_dir)
        print(f"Created project directory: {project_dir}")
    else:
        print(f"Using existing project directory: {project_dir}")
        
    # Create input directory
    input_dir = os.path.join(project_dir, "input")
    if not os.path.exists(input_dir):
        os.makedirs(input_dir)
        
    # Create logs directory
    logs_dir = os.path.join(project_dir, "logs")
    if not os.path.exists(logs_dir):
        os.makedirs(logs_dir)
        
    # Copy config file to input directory
    config_dest = os.path.join(input_dir, os.path.basename(args.config))
    if os.path.abspath(args.config) != os.path.abspath(config_dest):
        import shutil
        shutil.copy(args.config, config_dest)
        print(f"Copied config to {config_dest}")
    
    # Define log files
    log_step1 = os.path.join(logs_dir, "step1_pesto.log")
    log_step2 = os.path.join(logs_dir, "step2_haddock.log")
    log_step3 = os.path.join(logs_dir, "step3_boltz.log")
    log_step4 = os.path.join(logs_dir, "step4_double_pass.log")
    
    # Step 1
    if args.step == 0 or args.step == 1:
        target_pesto_pdb = step1_pesto(config, project_dir, input_dir, log_step1)
    else:
        # If skipping step 1, try to find the output
        step1_dir = os.path.join(project_dir, "step1_pesto")
        target_pdb = config['target_pdb']
        if os.path.exists(target_pdb):
            base = os.path.splitext(os.path.basename(target_pdb))[0]
        else:
            base = target_pdb
        target_pesto_pdb = os.path.join(step1_dir, f"{base}_pesto.pdb")
        
        if not os.path.exists(target_pesto_pdb):
             print(f"Warning: {target_pesto_pdb} does not exist. Step 2 might fail.")

    # Step 2
    if args.step == 0 or args.step == 2:
        haddock_run_dir = step2_haddock(config, target_pesto_pdb, project_dir, log_step2)
    else:
        step2_dir = os.path.join(project_dir, "step2_haddock")
        run_dir_name = config.get('haddock_run_dir', 'haddock_results')
        haddock_run_dir = os.path.join(step2_dir, os.path.basename(run_dir_name))

    # Step 3
    if args.step == 0 or args.step == 3:
        # Use the config in input_dir for step 3 to ensure consistency
        step3_boltz(config, haddock_run_dir, config_dest, project_dir, log_step3,
                    override=args.override, ref_pdb=abs_ref)

    # Step 4 (double-pass refinement)
    if args.step == 4 or (args.step == 0 and args.doublepass):
        step4_double_pass(config, config_dest, project_dir, log_step4,
                          override=args.override, ref_pdb=abs_ref)

if __name__ == "__main__":
    main()
