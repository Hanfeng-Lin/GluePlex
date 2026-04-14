import matplotlib
matplotlib.use("Agg")   # non-interactive backend (safe for subprocess calls)
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from sklearn.decomposition import PCA
import pandas as pd
import argparse

from bokeh.plotting import figure, output_file, save
from bokeh.models import ColumnDataSource, HoverTool, ColorBar, LinearColorMapper
from bokeh.transform import factor_cmap, linear_cmap
from bokeh.palettes import Spectral11, RdYlBu11

# ── Args ─────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="PCA visualization.")
parser.add_argument("--ligand_name",  default="Ligand",   help="Name of the ligand for plot titles.")
parser.add_argument("--protein_name", default="Protein",  help="Name of the protein for plot titles.")
args = parser.parse_args()

ligand_name  = args.ligand_name
protein_name = args.protein_name

# ── Load data ─────────────────────────────────────────────────────────────────
coordinates = pd.read_csv("atom.csv")
coordinates_data = coordinates.drop(
    columns=["ipTM", "file", "cluster", "pair_chains_ipTM_CRBN_VAV1_avg"]
).values

# ── PCA ───────────────────────────────────────────────────────────────────────
reducer  = PCA(n_components=2, random_state=42)
embedding = reducer.fit_transform(coordinates_data)

explained  = reducer.explained_variance_ratio_
pc1_label  = f"PC1 ({explained[0]*100:.1f}%)"
pc2_label  = f"PC2 ({explained[1]*100:.1f}%)"

# Shared data source columns
clusters     = coordinates["cluster"].astype(str).values
iptm_vals    = coordinates["ipTM"].values
chainp_vals  = coordinates["pair_chains_ipTM_CRBN_VAV1_avg"].values
file_vals    = coordinates["file"].values

# ══════════════════════════════════════════════════════════════════════════════
# BOKEH INTERACTIVE PLOTS
# ══════════════════════════════════════════════════════════════════════════════

TOOLS = "pan,wheel_zoom,box_zoom,reset,save"
POINT_SIZE = 7

# ── 1. Cluster plot ───────────────────────────────────────────────────────────
unique_clusters = sorted(set(clusters))
n_clusters = len(unique_clusters)
# Pick a palette large enough; tile Spectral11 if needed
palette = (Spectral11 * ((n_clusters // 11) + 1))[:n_clusters]

source_cluster = ColumnDataSource(dict(
    x=[float(v) for v in embedding[:, 0]],
    y=[float(v) for v in embedding[:, 1]],
    cluster=list(clusters),
    conformation=list(file_vals),
    iptm=[float(v) for v in iptm_vals],
    chainpair_iptm=[float(v) for v in chainp_vals],
))

p_cluster = figure(
    title=f"{protein_name} – {ligand_name}  |  Co-folding from Cluster",
    width=700, height=600,
    tools=TOOLS,
    x_axis_label=pc1_label, y_axis_label=pc2_label,
)
p_cluster.circle(
    "x", "y",
    source=source_cluster,
    size=POINT_SIZE,
    color=factor_cmap("cluster", palette=palette, factors=unique_clusters),
    legend_field="cluster",
    alpha=0.85,
)
p_cluster.add_tools(HoverTool(tooltips=[
    ("conformation", "@conformation"),
    ("cluster",      "@cluster"),
    ("ipTM",         "@iptm{0.000}"),
    ("chainpair_ipTM", "@chainpair_iptm{0.000}"),
]))
p_cluster.legend.title = "Cluster"
p_cluster.legend.click_policy = "hide"

output_file("PCA_cluster.html")
save(p_cluster)
print("Saved PCA_cluster.html")

# ── 2. ipTM plot ──────────────────────────────────────────────────────────────
source_iptm = ColumnDataSource(dict(
    x=[float(v) for v in embedding[:, 0]],
    y=[float(v) for v in embedding[:, 1]],
    cluster=list(clusters),
    conformation=list(file_vals),
    iptm=[float(v) for v in iptm_vals],
    chainpair_iptm=[float(v) for v in chainp_vals],
))

iptm_mapper = LinearColorMapper(palette=list(reversed(RdYlBu11)),
                                low=0.0, high=1.05)

p_iptm = figure(
    title=f"{protein_name} – {ligand_name}  |  ipTM Score",
    width=700, height=600,
    tools=TOOLS,
    x_axis_label=pc1_label, y_axis_label=pc2_label,
)
p_iptm.circle(
    "x", "y",
    source=source_iptm,
    size=POINT_SIZE,
    color={"field": "iptm", "transform": iptm_mapper},
    alpha=0.85,
)
p_iptm.add_tools(HoverTool(tooltips=[
    ("conformation", "@conformation"),
    ("cluster",      "@cluster"),
    ("ipTM",         "@iptm{0.000}"),
    ("chainpair_ipTM", "@chainpair_iptm{0.000}"),
]))
cbar_iptm = ColorBar(color_mapper=iptm_mapper, label_standoff=10, title="ipTM",
                     width=12, location=(0, 0))
p_iptm.add_layout(cbar_iptm, "right")

output_file("PCA_ipTM.html")
save(p_iptm)
print("Saved PCA_ipTM.html")

# ── 3. chainpair_ipTM plot ────────────────────────────────────────────────────
source_cp = ColumnDataSource(dict(
    x=[float(v) for v in embedding[:, 0]],
    y=[float(v) for v in embedding[:, 1]],
    cluster=list(clusters),
    conformation=list(file_vals),
    iptm=[float(v) for v in iptm_vals],
    chainpair_iptm=[float(v) for v in chainp_vals],
))

cp_mapper = LinearColorMapper(palette=list(reversed(RdYlBu11)),
                              low=0.0, high=1.05)

p_cp = figure(
    title=f"{protein_name} – {ligand_name}  |  Chain-pair ipTM Score",
    width=700, height=600,
    tools=TOOLS,
    x_axis_label=pc1_label, y_axis_label=pc2_label,
)
p_cp.circle(
    "x", "y",
    source=source_cp,
    size=POINT_SIZE,
    color={"field": "chainpair_iptm", "transform": cp_mapper},
    alpha=0.85,
)
p_cp.add_tools(HoverTool(tooltips=[
    ("conformation", "@conformation"),
    ("cluster",      "@cluster"),
    ("ipTM",         "@iptm{0.000}"),
    ("chainpair_ipTM", "@chainpair_iptm{0.000}"),
]))
cbar_cp = ColorBar(color_mapper=cp_mapper, label_standoff=10,
                   title="chainpair_ipTM", width=12, location=(0, 0))
p_cp.add_layout(cbar_cp, "right")

output_file("PCA_chainpair_ipTM.html")
save(p_cp)
print("Saved PCA_chainpair_ipTM.html")

# ══════════════════════════════════════════════════════════════════════════════
# STATIC PNG PLOTS (matplotlib)
# ══════════════════════════════════════════════════════════════════════════════

cmap_static = plt.cm.RdYlBu.reversed()

# ── 1. Cluster PNG ────────────────────────────────────────────────────────────
int_clusters = sorted(set(clusters))
cmap_cat     = matplotlib.colormaps.get_cmap("Spectral").resampled(len(int_clusters))
cluster_color = {c: cmap_cat(i) for i, c in enumerate(int_clusters)}

plt.figure(figsize=(6, 6))
for c in int_clusters:
    mask = clusters == c
    plt.scatter(embedding[mask, 0], embedding[mask, 1],
                color=cluster_color[c], s=20, label=str(c))
plt.legend(title="Cluster", bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
plt.xlabel(pc1_label, fontsize=14)
plt.ylabel(pc2_label, fontsize=14)
plt.title(f"{protein_name} – {ligand_name}  Co-folding from Cluster", pad=16, fontsize=16)
plt.tight_layout()
plt.savefig("PCA_cluster.png", dpi=300)
plt.close()

# ── 2. ipTM PNG ───────────────────────────────────────────────────────────────
plt.figure(figsize=(7, 6))
sc = plt.scatter(embedding[:, 0], embedding[:, 1],
                 c=iptm_vals, cmap=cmap_static, s=20, vmin=0, vmax=1.05)
cbar = plt.colorbar(sc, boundaries=np.linspace(0, 1.05, 256))
cbar.set_label("ipTM", rotation=270, labelpad=16)
cbar.set_ticks(np.arange(0, 1.01, 0.2))
plt.xlabel(pc1_label, fontsize=14)
plt.ylabel(pc2_label, fontsize=14)
plt.title(f"{protein_name} – {ligand_name}  Co-folding Score", pad=16, fontsize=16)
plt.tight_layout()
plt.savefig("PCA_ipTM.png", dpi=300)
plt.close()

# ── 3. chainpair_ipTM PNG ────────────────────────────────────────────────────
plt.figure(figsize=(7, 6))
sc = plt.scatter(embedding[:, 0], embedding[:, 1],
                 c=chainp_vals, cmap=cmap_static, s=20, vmin=0, vmax=1.05)
cbar = plt.colorbar(sc, boundaries=np.linspace(0, 1.05, 256))
cbar.set_label("chainpair_ipTM", rotation=270, labelpad=16)
cbar.set_ticks(np.arange(0, 1.01, 0.2))
plt.xlabel(pc1_label, fontsize=14)
plt.ylabel(pc2_label, fontsize=14)
plt.title(f"{protein_name} – {ligand_name}  Co-folding Score", pad=16, fontsize=16)
plt.tight_layout()
plt.savefig("PCA_chainpair_ipTM.png", dpi=300)
plt.close()

print("All PCA plots saved.")
