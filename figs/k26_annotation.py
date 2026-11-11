from gallery import *

# Annotation: marker-based (dot plot) vs reference-based (map onto an atlas; what the atlas lacks gets mislabelled)
fig = plt.figure(figsize=(9.8, 4.5))
gs = fig.add_gridspec(1, 2, width_ratios=[1.1, 1], wspace=0.25, left=0.08, right=0.99, top=0.8, bottom=0.14)
DARK = "#5E6670"; GREYF = "#B8BEC7"
from matplotlib.colors import LinearSegmentedColormap
CM = LinearSegmentedColormap.from_list("e", ["#E3E6EB", tint(0.45), CRIMSON])
genes = ["CD3E", "CD8A", "NKG7", "MS4A1", "LYZ", "COL1A1"]
clusters = ["1", "2", "3", "4", "5", "6"]
names = ["CD4 T", "CD8 T", "NK", "B", "monocyte", "fibroblast"]
pct = np.array([[90, 8, 10, 2, 3, 1], [88, 85, 40, 2, 2, 1], [10, 15, 92, 1, 4, 1], [3, 1, 2, 90, 5, 1],
                [2, 1, 6, 3, 95, 3], [1, 1, 1, 1, 6, 88]]) / 100
avg = pct * np.array([[1.0, 0.3, 0.3, 0.2, 0.2, 0.2], [1.0, 1.0, 0.5, 0.2, 0.2, 0.2], [0.3, 0.3, 1.0, 0.2, 0.3, 0.2],
                      [0.2, 0.2, 0.2, 1.0, 0.3, 0.2], [0.2, 0.2, 0.3, 0.2, 1.0, 0.3], [0.2, 0.2, 0.2, 0.2, 0.3, 1.0]])
a = fig.add_subplot(gs[0])
for i in range(6):
    for j in range(6):
        a.scatter(j, 5 - i, s=20 + 330 * pct[i, j], c=[avg[i, j]], cmap=CM, vmin=0, vmax=1, lw=0)
a.set_xticks(range(6)); a.set_xticklabels(genes, fontsize=8.4, color=TEXT, fontstyle="italic")
a.set_yticks(range(6)); a.set_yticklabels([f"cluster {c}" for c in clusters][::-1], fontsize=8.4, color=TEXT2)
for i, nm in enumerate(names):
    a.text(6.0, 5 - i, f"→ {nm}", ha="left", va="center", fontsize=9, color=TEXT, fontweight="semibold")
a.set_xlim(-0.6, 7.6); a.set_ylim(-0.6, 5.6)
for sp in ("left", "bottom"): a.spines[sp].set_visible(False)
a.tick_params(length=0)
a.set_title("(a) Marker-based", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
a.text(0, 1.03, "dot size = % of cells expressing; colour = mean level", transform=a.transAxes, fontsize=8.8, color=TEXT2)
a.text(0, -0.15, "✓ transparent   ✕ circular: you find the types you looked for", transform=a.transAxes, fontsize=8.6,
       color=TEXT)

b = fig.add_subplot(gs[1]); rng = np.random.default_rng(26)
ref = {"T cell": ([-1.6, 1.0], DARK), "B cell": ([1.4, 1.3], GREYF), "monocyte": ([0.0, -1.4], tint(0.45))}
for nm, (c, col) in ref.items():
    pts = rng.normal(c, 0.38, (120, 2)); b.scatter(pts[:, 0], pts[:, 1], s=6, color=col, lw=0, alpha=0.55)
    b.text(c[0], c[1] + 0.95, nm, ha="center", va="center", fontsize=8.6, color=TEXT2)
q = rng.normal([-1.45, 0.85], 0.22, (12, 2)); b.scatter(q[:, 0], q[:, 1], s=16, color=TEXT, lw=0, marker="^")
new = rng.normal([1.3, -1.0], 0.2, (12, 2)); b.scatter(new[:, 0], new[:, 1], s=18, color=CRIMSON, lw=0, marker="^")
b.annotate("", (0.35, -1.25), (1.05, -1.05), arrowprops=dict(arrowstyle="-|>", color=CRIMSON, lw=1, mutation_scale=9))
b.text(1.3, -1.55, "a type the atlas lacks →\nforced onto the nearest label", ha="center", va="top", fontsize=8.4,
       color=CRIMSON, fontweight="semibold", linespacing=1.15)
b.text(-0.95, 0.35, "query T cells:\nlabel transferred ✓", ha="left", va="top", fontsize=8.2, color=TEXT, linespacing=1.15)
b.set_xlim(-2.6, 2.6); b.set_ylim(-2.6, 2.3); b.set_xticks([]); b.set_yticks([])
for sp in b.spines.values(): sp.set_visible(True); sp.set_color("#E3E6EB")
b.set_title("(b) Reference-based", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
b.text(0, 1.03, "map your cells (▲) onto a labelled atlas (•)", transform=b.transAxes, fontsize=8.8, color=TEXT2)
b.text(0, -0.08, "✓ fast, consistent   ✕ inherits the atlas's names and gaps", transform=b.transAxes, fontsize=8.6,
       color=TEXT, va="top")
fig.text(0.08, -0.06, "Use both. Where they disagree, you have either a new population or an artefact.", fontsize=8.9, color=TEXT2)
header(a, "Single-cell RNA-seq · annotation", "Naming clusters: from markers, or from a reference", dy=40, x1=2.2)
fig.savefig(OUT + "annotation.png"); print("ok k26 annotation")
