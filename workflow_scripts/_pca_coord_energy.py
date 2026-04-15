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
from matplotlib.colors import rgb2hex

# Smooth continuous palette (256 interpolated colors) for ipTM / chainpair_ipTM.
# Bokeh's built-in RdYlBu11 is only 11 bins; interpolate via matplotlib for a
# truly smooth gradient in the HTML plots.
_cmap_continuous = plt.cm.RdYlBu.reversed()
CONTINUOUS_PALETTE = tuple(rgb2hex(_cmap_continuous(i / 255)) for i in range(256))

# ── Args ─────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="PCA visualization.")
parser.add_argument("--ligand_name",  default="Ligand",   help="Name of the ligand for plot titles.")
parser.add_argument("--protein_name", default="Protein",  help="Name of the protein for plot titles.")
parser.add_argument("--input_csv",    default="atom.csv", help="Input CSV written by 3.2.post_analysis.py.")
args = parser.parse_args()

ligand_name  = args.ligand_name
protein_name = args.protein_name

# ── Load data ─────────────────────────────────────────────────────────────────
coordinates = pd.read_csv(args.input_csv)
meta_cols = ["ipTM", "file", "cluster", "pair_chains_ipTM_avg"]
if "pass" in coordinates.columns:
    meta_cols.append("pass")
coordinates_data = coordinates.drop(columns=meta_cols).values

# ── PCA ───────────────────────────────────────────────────────────────────────
reducer  = PCA(n_components=2, random_state=42)
embedding = reducer.fit_transform(coordinates_data)

explained  = reducer.explained_variance_ratio_
pc1_label  = f"PC1 ({explained[0]*100:.1f}%)"
pc2_label  = f"PC2 ({explained[1]*100:.1f}%)"

# Shared data source columns
clusters     = coordinates["cluster"].astype(str).values
iptm_vals    = coordinates["ipTM"].values
chainp_vals  = coordinates["pair_chains_ipTM_avg"].values
file_vals    = coordinates["file"].values
if "pass" in coordinates.columns:
    pass_vals = coordinates["pass"].astype(int).values
else:
    pass_vals = np.ones(len(coordinates), dtype=int)

# Bokeh per-point edge styling based on pass (pass-2 gets a black ring)
edge_color = np.where(pass_vals == 2, "black", None)
edge_width = np.where(pass_vals == 2, 1.0, 0.0).astype(float)

# ══════════════════════════════════════════════════════════════════════════════
# BOKEH INTERACTIVE PLOTS
# ══════════════════════════════════════════════════════════════════════════════

TOOLS = "pan,wheel_zoom,box_zoom,reset,save"
POINT_SIZE = 7

def _shared_source():
    return ColumnDataSource(dict(
        x=[float(v) for v in embedding[:, 0]],
        y=[float(v) for v in embedding[:, 1]],
        cluster=list(clusters),
        conformation=list(file_vals),
        iptm=[float(v) for v in iptm_vals],
        chainpair_iptm=[float(v) for v in chainp_vals],
        pass_num=[int(v) for v in pass_vals],
        line_color=[c if c is not None else "black" for c in edge_color],  # invisible when line_width=0
        line_width=list(edge_width),
    ))

# ── 1. Cluster plot ───────────────────────────────────────────────────────────
unique_clusters = sorted(set(clusters))
n_clusters = len(unique_clusters)
palette = (Spectral11 * ((n_clusters // 11) + 1))[:n_clusters]

source_cluster = _shared_source()

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
    fill_color=factor_cmap("cluster", palette=palette, factors=unique_clusters),
    line_color="line_color",
    line_width="line_width",
    legend_field="cluster",
    alpha=0.85,
)
p_cluster.add_tools(HoverTool(tooltips=[
    ("conformation", "@conformation"),
    ("cluster",      "@cluster"),
    ("pass",         "@pass_num"),
    ("ipTM",         "@iptm{0.000}"),
    ("chainpair_ipTM", "@chainpair_iptm{0.000}"),
]))
p_cluster.legend.title = "Cluster"
p_cluster.legend.click_policy = "hide"

output_file("PCA_cluster.html")
save(p_cluster)
print("Saved PCA_cluster.html")

# ── 2. ipTM plot ──────────────────────────────────────────────────────────────
source_iptm = _shared_source()
iptm_mapper = LinearColorMapper(palette=CONTINUOUS_PALETTE, low=0.0, high=1.05)

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
    fill_color={"field": "iptm", "transform": iptm_mapper},
    line_color="line_color",
    line_width="line_width",
    alpha=0.85,
)
p_iptm.add_tools(HoverTool(tooltips=[
    ("conformation", "@conformation"),
    ("cluster",      "@cluster"),
    ("pass",         "@pass_num"),
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
source_cp = _shared_source()
cp_mapper = LinearColorMapper(palette=CONTINUOUS_PALETTE, low=0.0, high=1.05)

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
    fill_color={"field": "chainpair_iptm", "transform": cp_mapper},
    line_color="line_color",
    line_width="line_width",
    alpha=0.85,
)
p_cp.add_tools(HoverTool(tooltips=[
    ("conformation", "@conformation"),
    ("cluster",      "@cluster"),
    ("pass",         "@pass_num"),
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
# STATIC PNG PLOTS (matplotlib) — pass-2 points rendered with black edge
# ══════════════════════════════════════════════════════════════════════════════

cmap_static = plt.cm.RdYlBu.reversed()
p1_mask = pass_vals == 1
p2_mask = pass_vals == 2

# ── 1. Cluster PNG ────────────────────────────────────────────────────────────
int_clusters = sorted(set(clusters))
cmap_cat     = matplotlib.colormaps.get_cmap("Spectral").resampled(len(int_clusters))
cluster_color = {c: cmap_cat(i) for i, c in enumerate(int_clusters)}

plt.figure(figsize=(6, 6))
for c in int_clusters:
    mask = clusters == c
    m1 = mask & p1_mask
    m2 = mask & p2_mask
    if m1.any():
        plt.scatter(embedding[m1, 0], embedding[m1, 1],
                    color=cluster_color[c], s=20, label=str(c),
                    edgecolors="none")
    if m2.any():
        plt.scatter(embedding[m2, 0], embedding[m2, 1],
                    color=cluster_color[c], s=20,
                    label=None if m1.any() else str(c),
                    edgecolors="black", linewidths=0.8)
plt.legend(title="Cluster", bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
plt.xlabel(pc1_label, fontsize=14)
plt.ylabel(pc2_label, fontsize=14)
plt.title(f"{protein_name} – {ligand_name}  Co-folding from Cluster", pad=16, fontsize=16)
plt.tight_layout()
plt.savefig("PCA_cluster.png", dpi=300)
plt.close()

def _scatter_continuous(values, label, outfile, title_suffix):
    plt.figure(figsize=(7, 6))
    if p1_mask.any():
        plt.scatter(embedding[p1_mask, 0], embedding[p1_mask, 1],
                    c=values[p1_mask], cmap=cmap_static, s=20, vmin=0, vmax=1.05,
                    edgecolors="none")
    sc = None
    if p2_mask.any():
        sc = plt.scatter(embedding[p2_mask, 0], embedding[p2_mask, 1],
                         c=values[p2_mask], cmap=cmap_static, s=20, vmin=0, vmax=1.05,
                         edgecolors="black", linewidths=0.8)
    if sc is None:
        sc = plt.scatter(embedding[p1_mask, 0], embedding[p1_mask, 1],
                         c=values[p1_mask], cmap=cmap_static, s=20, vmin=0, vmax=1.05,
                         edgecolors="none")
    cbar = plt.colorbar(sc, boundaries=np.linspace(0, 1.05, 256))
    cbar.set_label(label, rotation=270, labelpad=16)
    cbar.set_ticks(np.arange(0, 1.01, 0.2))
    plt.xlabel(pc1_label, fontsize=14)
    plt.ylabel(pc2_label, fontsize=14)
    plt.title(f"{protein_name} – {ligand_name}  {title_suffix}", pad=16, fontsize=16)
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    plt.close()

# ── 2. ipTM PNG ───────────────────────────────────────────────────────────────
_scatter_continuous(iptm_vals, "ipTM", "PCA_ipTM.png", "Co-folding Score")

# ── 3. chainpair_ipTM PNG ────────────────────────────────────────────────────
_scatter_continuous(chainp_vals, "chainpair_ipTM", "PCA_chainpair_ipTM.png", "Co-folding Score")

print("All PCA plots saved.")
