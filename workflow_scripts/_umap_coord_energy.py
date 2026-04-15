import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import umap
import os
import pandas as pd
import umap.plot
import argparse
import sys

parser = argparse.ArgumentParser(description="UMAP visualization.")
parser.add_argument("--ligand_name", default="Ligand", help="Name of the ligand for plot titles.")
parser.add_argument("--protein_name", default="Protein", help="Name of the protein for plot titles.")
parser.add_argument("--input_csv", default="atom.csv", help="Input CSV written by 3.2.post_analysis.py.")
args = parser.parse_args()

ligand_name = args.ligand_name
protein_name = args.protein_name

coordinates = pd.read_csv(args.input_csv)

meta_cols = ["ipTM", "file", "cluster", "pair_chains_ipTM_avg"]
if "pass" in coordinates.columns:
    meta_cols.append("pass")
coordinates_data = coordinates.drop(columns=meta_cols).values

if "pass" in coordinates.columns:
    pass_vals = coordinates["pass"].astype(int).values
else:
    pass_vals = np.ones(len(coordinates), dtype=int)
p1_mask = pass_vals == 1
p2_mask = pass_vals == 2

reducer = umap.UMAP(n_neighbors=10, metric="euclidean", min_dist=1, random_state=42)
embedding = reducer.fit(coordinates_data)

hover_data = pd.DataFrame({'conformation': coordinates["file"].values,
                           'cluster': coordinates["cluster"].values,
                           'pass': pass_vals,
                           'ipTM': coordinates["ipTM"].values,
                           'chainpair_ipTM': coordinates["pair_chains_ipTM_avg"].values})

# Interactive cluster plot (umap.plot.interactive does not support per-point edges;
# we still emit it as a reference, but the PNG below distinguishes passes visually)
umap.plot.output_file("UMAP_cluster.html")
cluster_plot = umap.plot.interactive(embedding, labels=coordinates["cluster"], hover_data=hover_data,
                          color_key_cmap="Spectral", background="white", point_size=5)
umap.plot.show(cluster_plot)

embedding_xy = reducer.fit_transform(coordinates_data)

# ── Cluster PNG with pass-aware edges ──
plt.figure(figsize=(7, 6))
ax = plt.gca()
int_clusters = sorted(set(coordinates["cluster"].astype(str).values))
cmap_cat = matplotlib.colormaps.get_cmap("Spectral").resampled(len(int_clusters))
cluster_color = {c: cmap_cat(i) for i, c in enumerate(int_clusters)}
clusters_str = coordinates["cluster"].astype(str).values
for c in int_clusters:
    mask = clusters_str == c
    m1 = mask & p1_mask
    m2 = mask & p2_mask
    if m1.any():
        ax.scatter(embedding_xy[m1, 0], embedding_xy[m1, 1],
                   color=cluster_color[c], s=20, label=str(c), edgecolors="none")
    if m2.any():
        ax.scatter(embedding_xy[m2, 0], embedding_xy[m2, 1],
                   color=cluster_color[c], s=20,
                   label=None if m1.any() else str(c),
                   edgecolors="black", linewidths=0.8)
ax.legend(title="Cluster", bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
ax.set_title(f"{protein_name} - {ligand_name} Co-folding from Cluster", pad=16, fontsize=16)
ax.set_xlabel("UMAP_1", fontsize=14)
ax.set_ylabel("UMAP_2", fontsize=14)
ax.xaxis.set_major_locator(ticker.AutoLocator())
ax.yaxis.set_major_locator(ticker.AutoLocator())
plt.tight_layout()
plt.savefig("UMAP_cluster.png", dpi=300)
plt.close()

cmap = plt.cm.RdYlBu.reversed()

# Interactive ipTM plot
umap.plot.output_file("UMAP_ipTM.html")
ipTM_plot = umap.plot.interactive(embedding, values=coordinates["ipTM"].values, hover_data=hover_data,
                          cmap=cmap, background="white", point_size=5)
umap.plot.show(ipTM_plot)

def _scatter_continuous(values, label, outfile, title_suffix):
    plt.figure(figsize=(7, 6))
    if p1_mask.any():
        plt.scatter(embedding_xy[p1_mask, 0], embedding_xy[p1_mask, 1],
                    c=values[p1_mask], cmap=cmap, s=20, vmin=0, vmax=1.05,
                    edgecolors="none")
    sc = None
    if p2_mask.any():
        sc = plt.scatter(embedding_xy[p2_mask, 0], embedding_xy[p2_mask, 1],
                         c=values[p2_mask], cmap=cmap, s=20, vmin=0, vmax=1.05,
                         edgecolors="black", linewidths=0.8)
    if sc is None:
        sc = plt.scatter(embedding_xy[p1_mask, 0], embedding_xy[p1_mask, 1],
                         c=values[p1_mask], cmap=cmap, s=20, vmin=0, vmax=1.05,
                         edgecolors="none")
    cbar = plt.colorbar(sc, boundaries=np.linspace(0, 1.05, 256))
    cbar.set_label(label, rotation=270, labelpad=16)
    cbar.set_ticks(np.arange(0, 1.01, 0.2))
    plt.xlabel("UMAP_1", fontsize=14)
    plt.ylabel("UMAP_2", fontsize=14)
    plt.title(f"{protein_name} - {ligand_name} {title_suffix}", pad=16, fontsize=16)
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    plt.close()

_scatter_continuous(coordinates["ipTM"].values, "ipTM", "UMAP_ipTM.png", "Co-folding Score")

# Interactive chainpair_ipTM plot
umap.plot.output_file("UMAP_chainpair_ipTM.html")
chainpair_ipTM_plot = umap.plot.interactive(reducer, values=coordinates["pair_chains_ipTM_avg"].values, hover_data=hover_data,
                          cmap=cmap, background="white", point_size=5)
umap.plot.show(chainpair_ipTM_plot)

_scatter_continuous(coordinates["pair_chains_ipTM_avg"].values, "chainpair_ipTM",
                    "UMAP_chainpair_ipTM.png", "Co-folding Score")
