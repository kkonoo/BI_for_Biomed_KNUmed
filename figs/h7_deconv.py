from gallery import *
from matplotlib.colors import LinearSegmentedColormap

# A change in cell-type mix looks like differential expression in bulk RNA-seq
fig, axes = plt.subplots(1, 3, figsize=(9.8, 4.4), gridspec_kw={"wspace": 0.62, "width_ratios": [1.0, 0.8, 1.0]})
DARK = "#5E6670"; GREYF = "#B8BEC7"
CM = LinearSegmentedColormap.from_list("cr", ["#FFFFFF", tint(0.35), CRIMSON])
ct = ["Neuron", "Astrocyte", "Microglia", "Oligodendrocyte"]
prof = np.array([[9.2, 1.1, 0.4, 1.9], [1.0, 8.6, 0.6, 2.2], [0.3, 0.8, 9.0, 0.5], [1.5, 2.0, 0.7, 8.8]])
p_ctrl = np.array([0.42, 0.26, 0.07, 0.25])
p_case = np.array([0.31, 0.28, 0.16, 0.25])

def title(ax, t, sub):
    ax.set_title(t, loc="left", fontsize=12.5, fontweight="semibold", color=TEXT, pad=24)
    ax.annotate(sub, (0, 1), xycoords="axes fraction", xytext=(0, 6), textcoords="offset points",
                ha="left", va="bottom", fontsize=9.3, color=TEXT2)

# ① profiles
ax = axes[0]
ax.imshow(prof, cmap=CM, aspect="auto", vmin=0, vmax=10)
ax.set_xticks(range(4)); ax.set_xticklabels([f"gene {k}" for k in range(1, 5)], fontsize=8.5, color=TEXT2)
ax.set_yticks(range(4)); ax.set_yticklabels(ct, fontsize=9, color=TEXT)
ax.tick_params(length=0)
for s in ax.spines.values(): s.set_visible(False)
for i in range(4):
    for j in range(4):
        ax.text(j, i, f"{prof[i, j]:.1f}", ha="center", va="center", fontsize=9,
                color="white" if prof[i, j] > 5 else TEXT)
title(ax, "① Cell-type profiles", "each gene marks one cell type")

# ② proportions
ax = axes[1]
cols = [CRIMSON, tint(0.35), DARK, "#D5D9DF"]
bl = np.zeros(2)
for i, (c, col) in enumerate(zip(ct, cols)):
    v = np.array([p_ctrl[i], p_case[i]])
    ax.bar([0, 1], v, bottom=bl, color=col, width=0.62, edgecolor="white", lw=1.5)
    for x in (0, 1):
        ax.text(x, bl[x] + v[x] / 2, f"{v[x] * 100:.0f}%", ha="center", va="center", fontsize=8.5,
                color="white" if col in (CRIMSON, DARK) else TEXT)
    ax.text(1.38, bl[1] + v[1] / 2, c, ha="left", va="center", fontsize=8.5,
            color=TEXT if col != "#D5D9DF" else TEXT2)
    bl += v
ax.set_xticks([0, 1]); ax.set_xticklabels(["controls", "cases"], fontsize=9, color=TEXT)
ax.set_ylim(0, 1); ax.set_yticks([]); ax.set_xlim(-0.45, 1.35)
for s in ("left", "top", "right"): ax.spines[s].set_visible(False)
ax.spines["bottom"].set_color(LINE); ax.tick_params(length=0)
title(ax, "② The mix changes", "neuron loss, more microglia")

# ③ bulk result
ax = axes[2]
b_ctrl = p_ctrl @ prof; b_case = p_case @ prof
X = np.arange(4); w = 0.36
ax.bar(X - w / 2, b_ctrl, w, color=GREYF, edgecolor="white")
ax.bar(X + w / 2, b_case, w, color=CRIMSON, edgecolor="white")
for i, (a, b) in enumerate(zip(b_ctrl, b_case)):
    fc = b / a
    big = abs(np.log2(fc)) > 0.25
    ax.text(i, max(a, b) + 0.2, f"{fc:.2f}×", ha="center", va="bottom", fontsize=9,
            fontweight="bold" if big else "normal", color=CRIMSON if big else TEXT2)
ax.set_xticks(X); ax.set_xticklabels(["gene 1\nneuron", "gene 2", "gene 3\nmicroglia", "gene 4"], fontsize=8.5,
                                     color=TEXT2)
ax.set_ylim(0, max(b_ctrl.max(), b_case.max()) * 1.45); ax.set_yticks([])
for s in ("left", "top", "right"): ax.spines[s].set_visible(False)
ax.spines["bottom"].set_color(LINE); ax.tick_params(length=0)
ax.text(-0.5, ax.get_ylim()[1] * 0.97, "■", color=GREYF, fontsize=10, ha="left", va="top")
ax.text(-0.15, ax.get_ylim()[1] * 0.97, "controls", color=TEXT2, fontsize=8.8, ha="left", va="top")
ax.text(1.0, ax.get_ylim()[1] * 0.97, "■", color=CRIMSON, fontsize=10, ha="left", va="top")
ax.text(1.35, ax.get_ylim()[1] * 0.97, "cases", color=CRIMSON, fontsize=8.8, ha="left", va="top")
title(ax, "③ Bulk RNA-seq result", "'neuronal genes down, microglial up'")
fig.text(0.47, -0.06, "Not one cell changed its expression — only the proportions changed.", ha="center",
         va="center", fontsize=9.5, color=TEXT, fontweight="semibold")

header(axes[0], "Beyond DE · cell-type composition",
       "A change in the cell mix looks exactly like differential expression", x=-0.5, dy=44, x1=4.35)
fig.savefig(OUT + "deconvolution.png"); print("ok h7")
