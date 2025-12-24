# GluePlex

## Stepwise
1. PeSTO (TODO)
   
   https://pesto.epfl.ch/
   Download the generated TARGET_pesto.pdb to `2.HADDOCK3` folder
   
2. HADDOCK3
   Modify the filename in each script.
   Usage:
   ```
   python 0.fetch_active_res.py TARGET_pesto.pdb
   # Active residues with b-factor>=0.3 will show. Paste this string into 1.generate_actpass.py
   ```
   For CRBN-pesto.pdb, use the default active residues.

   ```
   python 1.generate_actpass.py
   # constraint.tbl will be generated
   ```
   ```
   python 2.run_haddock_workflow.py
   # Main HADDOCK3 docking workflow. Typical run 1 hr.
   ```
   ```
   python 3.export_mmcif.py
   # convert docking result into mmCIF using ChimeraX. All other converting software will fail in boltz2.
   ```
   Move all the `.cif` files to `3.Boltz2` folder.
   
3. Boltz2-templated co-folding
   Modify the filename in `1.run_batch_clusters.py`.
   ```
   python 1.run_batch_clusters.py  # Main boltz2 script
   python 2.post_analysis.py  # combine result and UMAP
   ```
