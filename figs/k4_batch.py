from gallery import *
rng = np.random.default_rng(4)

# Integration: colour by donor AND by cell type. Over-correction looks perfectly mixed by donor
# but merges a disease-only population into a neighbouring type.
fig, axes = plt.subplots(2, 3, figsize=(9.8, 6.4), gridspec_kw={"wspace": 0.06, "hspace": 0.08})
DON = ["#2B2F36", "#6B737D", "#A3AAB3", "#D3D7DC"]
TYP = ["#5E6670", "#A9B0B9", tint(0.4), CRIMSON]
centres = np.array([[-1.8, 1.2], [1.6, 1.4], [0.0, -1.5], [1.9, -1.4]])
cells = []
for d in range(4):
    for ty in range(4):
        if ty == 3 and d < 2: continue                 # type D only exists in the two patients
        n = 70 if ty < 3 else 55
        cells += [(d, ty, p) for p in rng.normal(0, 0.32, (n, 2))]
don = np.array([c[0] for c in cells]); typ = np.array([c[1] for c in cells]); noise = np.array([c[2] for c in cells])
shift = np.array([[-1.2, 1.0], [1.3, 0.8], [-1.0, -1.1], [1.1, -0.9]]) * 1.25
before = (centres[typ] * 0.55 + shift[don]) + noise * 0.8
good = centres[typ] + noise
over = good.copy(); over[typ == 3] = centres[2] + noise[typ == 3] * 1.05 + np.array([0.35, 0.1])
cols = [("Before integration", before, "cells group by donor, not biology"),
        ("Good integration", good, "donors mixed, disease population kept"),
        ("Over-corrected", over, "donors mixed, disease population erased")]
for j, (t, E, sub) in enumerate(cols):
    for i, (lab, key, pal) in enumerate((("by donor", don, DON), ("by cell type", typ, TYP))):
        ax = axes[i, j]
        o = rng.permutation(len(E))
        ax.scatter(E[o, 0], E[o, 1], s=6, c=np.array(pal)[key[o]], lw=0, alpha=0.9, rasterized=True)
        ax.set_xlim(-3.6, 3.6); ax.set_ylim(-3.1, 3.1); ax.set_xticks([]); ax.set_yticks([])
        bad = j == 2 and i == 1
        for sp in ax.spines.values():
            sp.set_visible(True); sp.set_color(CRIMSON if bad else "#E3E6EB"); sp.set_linewidth(1.2 if bad else 0.8)
        ax.text(0.03, 0.96, lab, transform=ax.transAxes, fontsize=8.6, color=TEXT2, va="top", zorder=5,
                bbox=dict(fc="white", ec="none", pad=1.5, alpha=0.9))
        if i == 0:
            ax.set_title(t, loc="left", fontsize=12.3, fontweight="semibold", color=CRIMSON if j == 2 else TEXT, pad=22)
            ax.text(0, 1.03, sub, transform=ax.transAxes, fontsize=8.8, color=TEXT2)
axes[1, 1].annotate("disease-only type", (1.9, -1.4), xytext=(0.5, -2.75), fontsize=8.3, color=CRIMSON,
                    ha="center", arrowprops=dict(arrowstyle="-", color=CRIMSON, lw=0.7))
axes[1, 2].annotate("merged into a\nneighbouring type", (0.4, -1.4), xytext=(2.2, -2.5), fontsize=8.3, color=CRIMSON,
                    ha="center", linespacing=1.1, arrowprops=dict(arrowstyle="-", color=CRIMSON, lw=0.7))
from matplotlib.lines import Line2D
h1 = [Line2D([], [], marker="o", ls="", color=c, ms=5) for c in DON]
h2 = [Line2D([], [], marker="o", ls="", color=c, ms=5) for c in TYP]
fig.legend(h1 + h2, ["control 1", "control 2", "patient 1", "patient 2", "type A", "type B", "type C",
                     "type D (patients only)"], loc="lower center", ncol=8, frameon=False, fontsize=8,
           bbox_to_anchor=(0.5, 0.0), labelcolor=TEXT, handletextpad=0.1, columnspacing=0.8)
fig.text(0.5, -0.03, "Both right-hand panels look equally well mixed by donor; only the cell-type colouring shows the loss",
         ha="center", fontsize=9, color=TEXT2)
header(axes[0, 0], "Single-cell RNA-seq · integration", "Always colour by both: donor and cell type", dy=40, x1=3.12)
fig.savefig(OUT + "integration_overcorrection.png"); print("ok k4 integration")
