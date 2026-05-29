import matplotlib
matplotlib.use("Agg")   # non-interactive backend (safe for subprocess calls)
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from sklearn.decomposition import PCA
import pandas as pd
import argparse

from bokeh.plotting import figure, output_file, save
from bokeh.models import ColumnDataSource, HoverTool, ColorBar, LinearColorMapper, LabelSet
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
parser.add_argument("--protein_name", default="Protein",  help="Name of the target protein (B) for plot titles.")
parser.add_argument("--protein_a_name", default="",       help="Name of the fixed protein (A, e.g. CRBN).")
parser.add_argument("--input_csv",    default="atom.csv", help="Input CSV written by 3.2.post_analysis.py.")
args = parser.parse_args()

ligand_name  = args.ligand_name
protein_name = args.protein_name
# Filename suffix: e.g. "_CRBN_G3BP2"
name_suffix = f"_{args.protein_a_name}_{protein_name}" if args.protein_a_name else ""

# ── Load data ─────────────────────────────────────────────────────────────────
coordinates = pd.read_csv(args.input_csv)
meta_cols = ["ipTM", "file", "cluster", "pair_chains_ipTM_avg"]
if "pass" in coordinates.columns:
    meta_cols.append("pass")
if "ref_rmsd" in coordinates.columns:
    meta_cols.append("ref_rmsd")
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

# Top-10 models by pair_chains_ipTM_avg (descending) for rank labels
top10_idx = np.argsort(-chainp_vals)[:10]

def _add_top10_labels(fig):
    """Annotate top-10 points with rank numbers 1..10."""
    label_source = ColumnDataSource(dict(
        x=[float(embedding[i, 0]) for i in top10_idx],
        y=[float(embedding[i, 1]) for i in top10_idx],
        rank=[str(r) for r in range(1, len(top10_idx) + 1)],
    ))
    fig.add_layout(LabelSet(
        x="x", y="y", text="rank", source=label_source,
        x_offset=6, y_offset=6, text_font_size="11pt",
        text_color="black", text_font_style="bold",
    ))

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
_add_top10_labels(p_cluster)

output_file(f"PCA_cluster{name_suffix}.html")
save(p_cluster)
print(f"Saved PCA_cluster{name_suffix}.html")

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
_add_top10_labels(p_iptm)

output_file(f"PCA_ipTM{name_suffix}.html")
save(p_iptm)
print(f"Saved PCA_ipTM{name_suffix}.html")

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
_add_top10_labels(p_cp)

output_file(f"PCA_chainpair_ipTM{name_suffix}.html")
save(p_cp)
print(f"Saved PCA_chainpair_ipTM{name_suffix}.html")

# ══════════════════════════════════════════════════════════════════════════════
# STATIC PNG PLOTS (matplotlib) — pass-2 points rendered with black edge
# ══════════════════════════════════════════════════════════════════════════════

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

cmap_static = plt.cm.RdYlBu.reversed()
p1_mask = pass_vals == 1
p2_mask = pass_vals == 2
MARKER_SIZE = 36  # publication-quality marker size

def _pass2_rank(c):
    """Extract pass-2 rank from cluster label 'p2_topN'; returns 0 if not pass-2."""
    if isinstance(c, str) and c.startswith("p2_top"):
        try:
            return int(c[6:])
        except ValueError:
            return 0
    return 0

# Plot order for pass-2 points: top10 → top1 so top1 stays on top
_p2_ranks = np.array([_pass2_rank(c) for c in clusters])
_p2_idx = np.where(p2_mask)[0]
p2_order = _p2_idx[np.argsort(-_p2_ranks[_p2_idx])] if len(_p2_idx) else _p2_idx

# ── 1. Cluster PNG ────────────────────────────────────────────────────────────
all_clusters = sorted(set(clusters))
p1_clusters = [c for c in all_clusters if not (isinstance(c, str) and c.startswith("p2_top"))]
p2_clusters = sorted(
    [c for c in all_clusters if isinstance(c, str) and c.startswith("p2_top")],
    key=lambda c: -_pass2_rank(c),  # top10 first, top1 last
)
int_clusters = p1_clusters + p2_clusters
cmap_cat     = matplotlib.colormaps.get_cmap("Spectral").resampled(len(int_clusters))
cluster_color = {c: cmap_cat(i) for i, c in enumerate(int_clusters)}

fig, ax = plt.subplots(figsize=(7, 6))
for c in int_clusters:
    mask = clusters == c
    m1 = mask & p1_mask
    m2 = mask & p2_mask
    if m1.any():
        ax.scatter(embedding[m1, 0], embedding[m1, 1],
                   color=cluster_color[c], s=MARKER_SIZE, label=str(c),
                   edgecolors="none")
    if m2.any():
        ax.scatter(embedding[m2, 0], embedding[m2, 1],
                   color=cluster_color[c], s=MARKER_SIZE,
                   label=None if m1.any() else str(c),
                   edgecolors="black", linewidths=0.8)
ax.legend(title="Cluster", bbox_to_anchor=(1.05, 1), loc="upper left",
          fontsize=10, title_fontsize=12)
ax.set_xlabel(pc1_label)
ax.set_ylabel(pc2_label)
ax.set_title(f"{protein_name} – {ligand_name}  Co-folding from Cluster", pad=12)
fig.tight_layout()
fig.savefig(f"PCA_cluster{name_suffix}.png", dpi=600, bbox_inches="tight")
plt.close(fig)

def _scatter_continuous(values, label, outfile, title_suffix):
    fig, ax = plt.subplots(figsize=(7, 6))
    if p1_mask.any():
        ax.scatter(embedding[p1_mask, 0], embedding[p1_mask, 1],
                   c=values[p1_mask], cmap=cmap_static, s=MARKER_SIZE, vmin=0, vmax=1.05,
                   edgecolors="none")
    sc = None
    if len(p2_order):
        sc = ax.scatter(embedding[p2_order, 0], embedding[p2_order, 1],
                        c=values[p2_order], cmap=cmap_static, s=MARKER_SIZE, vmin=0, vmax=1.05,
                        edgecolors="black", linewidths=0.8)
    if sc is None:
        sc = ax.scatter(embedding[p1_mask, 0], embedding[p1_mask, 1],
                        c=values[p1_mask], cmap=cmap_static, s=MARKER_SIZE, vmin=0, vmax=1.05,
                        edgecolors="none")
    cbar = fig.colorbar(sc, ax=ax, boundaries=np.linspace(0, 1.05, 256))
    cbar.set_label(label, rotation=270, labelpad=18, fontsize=14)
    cbar.set_ticks(np.arange(0, 1.01, 0.2))
    cbar.ax.tick_params(labelsize=12)
    ax.set_xlabel(pc1_label)
    ax.set_ylabel(pc2_label)
    ax.set_title(f"{protein_name} – {ligand_name}  {title_suffix}", pad=12)
    fig.tight_layout()
    fig.savefig(outfile, dpi=600, bbox_inches="tight")
    plt.close(fig)

# ── 2. ipTM PNG ───────────────────────────────────────────────────────────────
_scatter_continuous(iptm_vals, "ipTM", f"PCA_ipTM{name_suffix}.png", "Co-folding Score")

# ── 3. chainpair_ipTM PNG ────────────────────────────────────────────────────
_scatter_continuous(chainp_vals, "chainpair_ipTM", f"PCA_chainpair_ipTM{name_suffix}.png", "Co-folding Score")

# ══════════════════════════════════════════════════════════════════════════════
# REF RMSD PLOTS (only if ref_rmsd column exists and has values)
# ══════════════════════════════════════════════════════════════════════════════

if "ref_rmsd" in coordinates.columns and coordinates["ref_rmsd"].notna().any():
    rmsd_vals = coordinates["ref_rmsd"].fillna(0).values.astype(float)
    rmsd_clamped = np.clip(rmsd_vals, 0, 10)

    # Bokeh interactive
    _cmap_rmsd_bokeh = plt.cm.RdYlBu
    _rmsd_palette = tuple(rgb2hex(_cmap_rmsd_bokeh(i / 255)) for i in range(256))
    rmsd_mapper = LinearColorMapper(palette=_rmsd_palette, low=0.0, high=10.0)
    source_rmsd = _shared_source()
    source_rmsd.data["ref_rmsd"] = [float(v) for v in rmsd_clamped]

    p_rmsd = figure(
        title=f"{protein_name} – {ligand_name}  |  Ref RMSD (chain B CA)",
        width=700, height=600, tools=TOOLS,
        x_axis_label=pc1_label, y_axis_label=pc2_label,
    )
    p_rmsd.circle(
        "x", "y", source=source_rmsd, size=POINT_SIZE,
        fill_color={"field": "ref_rmsd", "transform": rmsd_mapper},
        line_color="line_color", line_width="line_width", alpha=0.85,
    )
    p_rmsd.add_tools(HoverTool(tooltips=[
        ("conformation", "@conformation"), ("cluster", "@cluster"),
        ("pass", "@pass_num"), ("ref_rmsd", "@ref_rmsd{0.00}"),
        ("ipTM", "@iptm{0.000}"), ("chainpair_ipTM", "@chainpair_iptm{0.000}"),
    ]))
    cbar_rmsd = ColorBar(color_mapper=rmsd_mapper, label_standoff=10,
                         title="ref_rmsd (Å)", width=12, location=(0, 0))
    p_rmsd.add_layout(cbar_rmsd, "right")
    _add_top10_labels(p_rmsd)
    output_file(f"PCA_ref_rmsd{name_suffix}.html")
    save(p_rmsd)
    print(f"Saved PCA_ref_rmsd{name_suffix}.html")

    # Static PNG (0-10 range)
    cmap_rmsd = plt.cm.RdYlBu
    fig, ax = plt.subplots(figsize=(7, 6))
    if p1_mask.any():
        ax.scatter(embedding[p1_mask, 0], embedding[p1_mask, 1],
                   c=rmsd_clamped[p1_mask], cmap=cmap_rmsd, s=MARKER_SIZE, vmin=0, vmax=10,
                   edgecolors="none")
    sc = None
    if len(p2_order):
        sc = ax.scatter(embedding[p2_order, 0], embedding[p2_order, 1],
                        c=rmsd_clamped[p2_order], cmap=cmap_rmsd, s=MARKER_SIZE, vmin=0, vmax=10,
                        edgecolors="black", linewidths=0.8)
    if sc is None:
        sc = ax.scatter(embedding[p1_mask, 0], embedding[p1_mask, 1],
                        c=rmsd_clamped[p1_mask], cmap=cmap_rmsd, s=MARKER_SIZE, vmin=0, vmax=10,
                        edgecolors="none")
    cbar = fig.colorbar(sc, ax=ax, boundaries=np.linspace(0, 10, 256))
    cbar.set_label("ref_rmsd (Å)", rotation=270, labelpad=18, fontsize=14)
    cbar.set_ticks(np.arange(0, 10.1, 2))
    cbar.ax.tick_params(labelsize=12)
    ax.set_xlabel(pc1_label)
    ax.set_ylabel(pc2_label)
    ax.set_title(f"{protein_name} – {ligand_name}  Ref RMSD", pad=12)
    fig.tight_layout()
    fig.savefig(f"PCA_ref_rmsd{name_suffix}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved PCA_ref_rmsd{name_suffix}.png")

print("All PCA plots saved.")
