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
from bokeh.models import ColumnDataSource, LabelSet

parser = argparse.ArgumentParser(description="UMAP visualization.")
parser.add_argument("--ligand_name", default="Ligand", help="Name of the ligand for plot titles.")
parser.add_argument("--protein_name", default="Protein", help="Name of the target protein (B) for plot titles.")
parser.add_argument("--protein_a_name", default="", help="Name of the fixed protein (A, e.g. CRBN).")
parser.add_argument("--input_csv", default="atom.csv", help="Input CSV written by 3.2.post_analysis.py.")
args = parser.parse_args()

ligand_name = args.ligand_name
protein_name = args.protein_name
name_suffix = f"_{args.protein_a_name}_{protein_name}" if args.protein_a_name else ""

coordinates = pd.read_csv(args.input_csv)

meta_cols = ["ipTM", "file", "cluster", "pair_chains_ipTM_avg"]
if "pass" in coordinates.columns:
    meta_cols.append("pass")
if "ref_rmsd" in coordinates.columns:
    meta_cols.append("ref_rmsd")
coordinates_data = coordinates.drop(columns=meta_cols).values

if "pass" in coordinates.columns:
    pass_vals = coordinates["pass"].astype(int).values
else:
    pass_vals = np.ones(len(coordinates), dtype=int)
p1_mask = pass_vals == 1
p2_mask = pass_vals == 2

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 12,
    "axes.labelsize": 16,
    "axes.titlesize": 18,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
    "axes.linewidth": 1.4,
    "xtick.major.width": 1.2,
    "ytick.major.width": 1.2,
    "xtick.major.size": 5,
    "ytick.major.size": 5,
})
MARKER_SIZE = 36

def _pass2_rank(c):
    """Extract pass-2 rank from cluster label 'p2_topN'; returns 0 if not pass-2."""
    if isinstance(c, str) and c.startswith("p2_top"):
        try:
            return int(c[6:])
        except ValueError:
            return 0
    return 0

# Plot order for pass-2 points: top10 → top1 so top1 stays on top
_clusters_str_all = coordinates["cluster"].astype(str).values
_p2_ranks = np.array([_pass2_rank(c) for c in _clusters_str_all])
_p2_idx = np.where(p2_mask)[0]
p2_order = _p2_idx[np.argsort(-_p2_ranks[_p2_idx])] if len(_p2_idx) else _p2_idx

reducer = umap.UMAP(n_neighbors=10, metric="euclidean", min_dist=1, random_state=42)
embedding = reducer.fit(coordinates_data)
embedding_xy = reducer.fit_transform(coordinates_data)

# Top-10 models by pair_chains_ipTM_avg (descending) for rank labels
_chainp_vals = coordinates["pair_chains_ipTM_avg"].values
top10_idx = np.argsort(-_chainp_vals)[:10]

def _add_top10_labels(fig):
    """Annotate top-10 points with rank numbers 1..10 on a Bokeh figure."""
    label_source = ColumnDataSource(dict(
        x=[float(embedding_xy[i, 0]) for i in top10_idx],
        y=[float(embedding_xy[i, 1]) for i in top10_idx],
        rank=[str(r) for r in range(1, len(top10_idx) + 1)],
    ))
    fig.add_layout(LabelSet(
        x="x", y="y", text="rank", source=label_source,
        x_offset=6, y_offset=6, text_font_size="11pt",
        text_color="black", text_font_style="bold",
    ))

hover_data = pd.DataFrame({'conformation': coordinates["file"].values,
                           'cluster': coordinates["cluster"].values,
                           'pass': pass_vals,
                           'ipTM': coordinates["ipTM"].values,
                           'chainpair_ipTM': coordinates["pair_chains_ipTM_avg"].values})

# Interactive cluster plot (umap.plot.interactive does not support per-point edges;
# we still emit it as a reference, but the PNG below distinguishes passes visually)
umap.plot.output_file(f"UMAP_cluster{name_suffix}.html")
cluster_plot = umap.plot.interactive(embedding, labels=coordinates["cluster"], hover_data=hover_data,
                          color_key_cmap="Spectral", background="white", point_size=5)
_add_top10_labels(cluster_plot)
umap.plot.show(cluster_plot)

# ── Cluster PNG with pass-aware edges ──
fig, ax = plt.subplots(figsize=(7, 6))
all_clusters = sorted(set(coordinates["cluster"].astype(str).values))
p1_clusters = [c for c in all_clusters if not c.startswith("p2_top")]
p2_clusters = sorted([c for c in all_clusters if c.startswith("p2_top")],
                     key=lambda c: -_pass2_rank(c))  # top10 first, top1 last
int_clusters = p1_clusters + p2_clusters
cmap_cat = matplotlib.colormaps.get_cmap("Spectral").resampled(len(int_clusters))
cluster_color = {c: cmap_cat(i) for i, c in enumerate(int_clusters)}
clusters_str = coordinates["cluster"].astype(str).values
for c in int_clusters:
    mask = clusters_str == c
    m1 = mask & p1_mask
    m2 = mask & p2_mask
    if m1.any():
        ax.scatter(embedding_xy[m1, 0], embedding_xy[m1, 1],
                   color=cluster_color[c], s=MARKER_SIZE, label=str(c), edgecolors="none")
    if m2.any():
        ax.scatter(embedding_xy[m2, 0], embedding_xy[m2, 1],
                   color=cluster_color[c], s=MARKER_SIZE,
                   label=None if m1.any() else str(c),
                   edgecolors="black", linewidths=0.8)
ax.legend(title="Cluster", bbox_to_anchor=(1.05, 1), loc="upper left",
          fontsize=10, title_fontsize=12)
ax.set_title(f"{protein_name} – {ligand_name}  Co-folding from Cluster", pad=12)
ax.set_xlabel("UMAP_1")
ax.set_ylabel("UMAP_2")
ax.xaxis.set_major_locator(ticker.AutoLocator())
ax.yaxis.set_major_locator(ticker.AutoLocator())
fig.tight_layout()
fig.savefig(f"UMAP_cluster{name_suffix}.png", dpi=600, bbox_inches="tight")
plt.close(fig)

cmap = plt.cm.RdYlBu.reversed()

# Interactive ipTM plot
umap.plot.output_file(f"UMAP_ipTM{name_suffix}.html")
ipTM_plot = umap.plot.interactive(embedding, values=coordinates["ipTM"].values, hover_data=hover_data,
                          cmap=cmap, background="white", point_size=5)
_add_top10_labels(ipTM_plot)
umap.plot.show(ipTM_plot)

def _scatter_continuous(values, label, outfile, title_suffix):
    fig, ax = plt.subplots(figsize=(7, 6))
    if p1_mask.any():
        ax.scatter(embedding_xy[p1_mask, 0], embedding_xy[p1_mask, 1],
                   c=values[p1_mask], cmap=cmap, s=MARKER_SIZE, vmin=0, vmax=1.05,
                   edgecolors="none")
    sc = None
    if len(p2_order):
        sc = ax.scatter(embedding_xy[p2_order, 0], embedding_xy[p2_order, 1],
                        c=values[p2_order], cmap=cmap, s=MARKER_SIZE, vmin=0, vmax=1.05,
                        edgecolors="black", linewidths=0.8)
    if sc is None:
        sc = ax.scatter(embedding_xy[p1_mask, 0], embedding_xy[p1_mask, 1],
                        c=values[p1_mask], cmap=cmap, s=MARKER_SIZE, vmin=0, vmax=1.05,
                        edgecolors="none")
    cbar = fig.colorbar(sc, ax=ax, boundaries=np.linspace(0, 1.05, 256))
    cbar.set_label(label, rotation=270, labelpad=18, fontsize=14)
    cbar.set_ticks(np.arange(0, 1.01, 0.2))
    cbar.ax.tick_params(labelsize=12)
    ax.set_xlabel("UMAP_1")
    ax.set_ylabel("UMAP_2")
    ax.set_title(f"{protein_name} – {ligand_name}  {title_suffix}", pad=12)
    fig.tight_layout()
    fig.savefig(outfile, dpi=600, bbox_inches="tight")
    plt.close(fig)

_scatter_continuous(coordinates["ipTM"].values, "ipTM", f"UMAP_ipTM{name_suffix}.png", "Co-folding Score")

# Interactive chainpair_ipTM plot
umap.plot.output_file(f"UMAP_chainpair_ipTM{name_suffix}.html")
chainpair_ipTM_plot = umap.plot.interactive(reducer, values=coordinates["pair_chains_ipTM_avg"].values, hover_data=hover_data,
                          cmap=cmap, background="white", point_size=5)
_add_top10_labels(chainpair_ipTM_plot)
umap.plot.show(chainpair_ipTM_plot)

_scatter_continuous(coordinates["pair_chains_ipTM_avg"].values, "chainpair_ipTM",
                    f"UMAP_chainpair_ipTM{name_suffix}.png", "Co-folding Score")

# ── Ref RMSD plots (if column exists and has values) ──
if "ref_rmsd" in coordinates.columns and coordinates["ref_rmsd"].notna().any():
    rmsd_vals = coordinates["ref_rmsd"].fillna(0).values.astype(float)
    rmsd_clamped = np.clip(rmsd_vals, 0, 10)

    hover_rmsd = pd.DataFrame({
        'conformation': coordinates["file"].values,
        'cluster': coordinates["cluster"].values,
        'pass': pass_vals,
        'ref_rmsd': rmsd_clamped,
        'ipTM': coordinates["ipTM"].values,
        'chainpair_ipTM': coordinates["pair_chains_ipTM_avg"].values,
    })
    umap.plot.output_file(f"UMAP_ref_rmsd{name_suffix}.html")
    cmap_rmsd_interactive = plt.cm.RdYlBu
    rmsd_plot = umap.plot.interactive(reducer, values=rmsd_clamped, hover_data=hover_rmsd,
                                      cmap=cmap_rmsd_interactive, background="white", point_size=5)
    _add_top10_labels(rmsd_plot)
    umap.plot.show(rmsd_plot)

    # Static PNG (0-10 clamped)
    cmap_rmsd = plt.cm.RdYlBu
    fig, ax = plt.subplots(figsize=(7, 6))
    if p1_mask.any():
        ax.scatter(embedding_xy[p1_mask, 0], embedding_xy[p1_mask, 1],
                   c=rmsd_clamped[p1_mask], cmap=cmap_rmsd, s=MARKER_SIZE, vmin=0, vmax=10,
                   edgecolors="none")
    sc = None
    if len(p2_order):
        sc = ax.scatter(embedding_xy[p2_order, 0], embedding_xy[p2_order, 1],
                        c=rmsd_clamped[p2_order], cmap=cmap_rmsd, s=MARKER_SIZE, vmin=0, vmax=10,
                        edgecolors="black", linewidths=0.8)
    if sc is None:
        sc = ax.scatter(embedding_xy[p1_mask, 0], embedding_xy[p1_mask, 1],
                        c=rmsd_clamped[p1_mask], cmap=cmap_rmsd, s=MARKER_SIZE, vmin=0, vmax=10,
                        edgecolors="none")
    cbar_r = fig.colorbar(sc, ax=ax, boundaries=np.linspace(0, 10, 256))
    cbar_r.set_label("ref_rmsd (Å)", rotation=270, labelpad=18, fontsize=14)
    cbar_r.set_ticks(np.arange(0, 10.1, 2))
    cbar_r.ax.tick_params(labelsize=12)
    ax.set_xlabel("UMAP_1")
    ax.set_ylabel("UMAP_2")
    ax.set_title(f"{protein_name} – {ligand_name}  Ref RMSD", pad=12)
    fig.tight_layout()
    fig.savefig(f"UMAP_ref_rmsd{name_suffix}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved UMAP_ref_rmsd{name_suffix}.png")
