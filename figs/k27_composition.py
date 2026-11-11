from gallery import *

# Compositional data: one cell type expands, and every other type's proportion falls with no change in its own number
fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.2), gridspec_kw={"wspace": 0.35})
DARK = "#5E6670"; GREYF = "#B8BEC7"
types = [("A", CRIMSON), ("B", DARK), ("C", GREYF)]
ctrl = np.array([400, 300, 300]); dis = np.array([1200, 300, 300])
for ax, data, title, sub, ylab in ((axes[0], (ctrl, dis), "(a) Cell numbers", "what happened in the tissue", "cells"),
                                   (axes[1], (ctrl / ctrl.sum(), dis / dis.sum()), "(b) Proportions",
                                    "what a per-sample fraction shows", "share of cells")):
    for x, v in enumerate(data):
        bottom = 0
        for (nm, c), val in zip(types, v):
            ax.bar(x, val, bottom=bottom, color=c, width=0.55, edgecolor="white", lw=1.2)
            lab = f"{nm} {val:.0%}" if ylab.startswith("share") else f"{nm} {val:,}"
            ax.text(x, bottom + val / 2, lab, ha="center", va="center", fontsize=8.8,
                    color="white" if c != GREYF else TEXT, fontweight="semibold")
            bottom += val
    ax.set_xticks([0, 1]); ax.set_xticklabels(["control", "disease"], fontsize=9.3, color=TEXT)
    ax.set_ylabel(ylab, color=TEXT2, fontsize=9)
    for sp in ("left", "bottom"): ax.spines[sp].set_color(LINE)
    ax.tick_params(colors=TEXT2, labelsize=8.3); ax.tick_params(axis="x", length=0)
    ax.set_title(title, loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
    ax.text(0, 1.03, sub, transform=ax.transAxes, fontsize=9, color=TEXT2)
axes[1].set_yticks([0, 0.5, 1]); axes[1].set_yticklabels(["0", "50%", "100%"])
axes[0].text(1.38, 1500, "only A changed\n(×3)", ha="left", va="center", fontsize=9, color=CRIMSON, fontweight="semibold",
             linespacing=1.2)
axes[1].text(1.38, 0.83, "B and C 'fall'\n30% → 17%\nwithout changing", ha="left", va="center", fontsize=9, color=DARK,
             fontweight="semibold", linespacing=1.2)
fig.text(0.125, -0.06, "Proportions are tied together: use compositional methods (scCODA, miloR, propeller), not a χ² test on "
         "cell counts.  n is still the number of donors.", fontsize=8.8, color=TEXT2)
header(axes[0], "Single-cell RNA-seq · differential abundance", "When one type grows, every other share shrinks",
       dy=40, x1=2.45)
fig.savefig(OUT + "composition.png"); print("ok k27 composition")
