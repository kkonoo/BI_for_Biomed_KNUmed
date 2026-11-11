from gallery import *
rng = np.random.default_rng(2)

# QC: each threshold that removes an artefact also removes a real population
fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.6), gridspec_kw={"wspace": 0.22})
DARK = "#5E6670"; GREYF = "#C3C8CF"
def pop(n, umi, sd, gexp, mito, msd):
    u = np.exp(rng.normal(np.log(umi), sd, n)); g = 0.9 * u ** gexp * np.exp(rng.normal(0, 0.08, n))
    m = np.clip(rng.normal(mito, msd, n), 0.2, 60); return u, g, m
# (data, colour in panel a, colour in panel b): red = a real population that this panel's cutoff removes
P = {"ordinary cells": (pop(1400, 5000, 0.45, 0.82, 2.8, 1.1), GREYF, GREYF),
     "doublets": (pop(70, 15000, 0.2, 0.82, 2.8, 1.0), DARK, DARK),
     "large, active cells": (pop(45, 16000, 0.18, 0.82, 3.0, 1.0), CRIMSON, GREYF),
     "debris / broken cells": (pop(170, 700, 0.45, 0.85, 21, 6), DARK, DARK),
     "platelets, small cells": (pop(90, 650, 0.25, 0.84, 2.5, 1.0), CRIMSON, GREYF),
     "mito-rich cells": (pop(110, 5000, 0.4, 0.82, 14, 3.5), GREYF, CRIMSON)}
a, b = axes
for nm, ((u, g, m), ca, cb) in P.items():
    for ax_, xx, yy, c in ((a, u, g, ca), (b, g, m, cb)):
        ax_.scatter(xx, yy, s=15 if c != GREYF else 9, color=c, lw=0, alpha=0.85 if c != GREYF else 0.7,
                    zorder=3 if c == CRIMSON else 2)
for x in (a, b):
    for sp in ("left", "bottom"): x.spines[sp].set_color(LINE)
    x.tick_params(colors=TEXT2, labelsize=8.3); x.minorticks_off()
a.set_xscale("log"); a.set_yscale("log"); b.set_xscale("log")
a.set_xlabel("total UMI per barcode", color=TEXT2, fontsize=9.3)
a.set_ylabel("genes detected", color=TEXT2, fontsize=9.3)
b.set_xlabel("genes detected", color=TEXT2, fontsize=9.3); b.set_ylabel("mitochondrial reads (%)", color=TEXT2, fontsize=9.3)
a.set_xlim(150, 4e4); a.set_ylim(90, 6e3); b.set_xlim(90, 6e3); b.set_ylim(0, 45)
lo, hi = 500, 2000
for yv in (lo, hi): a.axhline(yv, color=CRIMSON, lw=0.9, ls=(0, (4, 3)))
a.text(160, lo * 0.86, "min genes = 500", ha="left", va="top", fontsize=8.2, color=CRIMSON)
a.text(160, hi * 1.08, "max genes", ha="left", va="bottom", fontsize=8.2, color=CRIMSON)
b.axhline(5, color=CRIMSON, lw=0.9, ls=(0, (4, 3)))
b.text(95, 5.6, "mito cutoff 5%", ha="left", va="bottom", fontsize=8.2, color=CRIMSON)
def lab(ax, x, y, t, c, ha="left"):
    ax.text(x, y, t, ha=ha, va="center", fontsize=8.5, color=c, fontweight="semibold", linespacing=1.15)
lab(a, 3.6e4, 4300, "doublets\n+ large real cells", TEXT, "right")
lab(a, 1300, 190, "← debris + platelets,\n    small cells", TEXT, "left")
lab(b, 600, 38, "broken cells:\nhigh mito, few genes", DARK)
lab(b, 2300, 24, "mito-rich cells\n(cardiomyocytes)", CRIMSON, "center")
a.set_title("(a) Gene-count cutoffs", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
a.text(0, 1.03, "each line removes an artefact and a real cell type", transform=a.transAxes, fontsize=9, color=TEXT2)
b.set_title("(b) Mito % cutoff", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
b.text(0, 1.03, "read it together with the gene count", transform=b.transAxes, fontsize=9, color=TEXT2)
# legend
from matplotlib.lines import Line2D
hs = [Line2D([], [], marker="o", ls="", color=c, ms=6) for c in (GREYF, DARK, CRIMSON)]
fig.legend(hs, ["ordinary cells", "artefacts (should go)", "real cells the cutoff also removes"], loc="lower center",
           ncol=3, frameon=False, fontsize=8.8, bbox_to_anchor=(0.5, -0.07), labelcolor=TEXT, handletextpad=0.3)
header(a, "Single-cell RNA-seq · QC", "Every QC threshold also deletes a real population", dy=40, x1=2.22)
fig.savefig(OUT + "qc_thresholds.png"); print("ok k2 qc")
