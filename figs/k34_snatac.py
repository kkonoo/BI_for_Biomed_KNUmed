from gallery import *
rng = np.random.default_rng(7)

# snATAC-seq: (a) the cells × peaks matrix is very sparse and nearly binary; (b) a gene activity score adds up the
# accessibility over a gene's body and promoter in each cell, giving a gene-like value to match RNA markers against
fig, axes = plt.subplots(1, 2, figsize=(9.8, 3.75), gridspec_kw={"wspace": 0.1, "width_ratios": [1, 1.12]})
DARK = "#5E6670"; GREYF = "#E3E6EB"
for ax in axes:
    ax.set_xlim(0, 10); ax.set_ylim(0.25, 4.95); ax.axis("off")


def cell(ax, x, y, v, s=0.36):
    ax.add_patch(mp.Rectangle((x, y), s, s, fc=tint(v) if v > 0 else GREYF, ec="white", lw=1.2))


# (a) cells × peaks, mostly 0
a = axes[0]; nc, npk, s = 8, 11, 0.36; x0, y0 = 1.1, 1.4
M = (rng.random((nc, npk)) < 0.13).astype(float)
M[rng.random((nc, npk)) < 0.03] = 2
for i in range(nc):
    for j in range(npk):
        cell(a, x0 + j * s, y0 + (nc - 1 - i) * s, {0: 0, 1: 0.55, 2: 1.0}[M[i, j]])
a.text(x0 + npk * s / 2, y0 + nc * s + 0.25, "peaks (open regions)", ha="center", fontsize=8.8, color=TEXT2)
a.text(x0 - 0.15, y0 + nc * s / 2, "cells", ha="right", va="center", fontsize=8.8, color=TEXT2, rotation=90)
lx = x0 + npk * s + 0.35
for k, (v, lab) in enumerate(((0, "0"), (0.55, "1"), (1.0, "2"))):
    cell(a, lx, y0 + 2.2 - k * 0.5, v)
    a.text(lx + 0.5, y0 + 2.2 - k * 0.5 + s / 2, lab + (" read" if lab == "1" else " reads"), va="center",
           fontsize=8.6, color=TEXT2)
a.text(0.4, 0.75, "Each locus has only two copies per cell → counts are nearly binary.\n"
       "Analysed as open vs. closed: LSI instead of PCA (ArchR, Signac)", ha="left", va="center",
       fontsize=8.6, color=TEXT, linespacing=1.35)
a.set_title("(a) A very sparse, nearly binary matrix", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT,
            pad=24)
a.text(0, 1.03, f"{100 * (M == 0).mean():.0f}% zeros here; ~95–99% in real data", transform=a.transAxes, ha="left",
       va="bottom", fontsize=9.1, color=TEXT2)

# (b) gene activity score
b = axes[1]; gx0, gx1, tss = 3.0, 8.6, 3.6                      # gene body, TSS on the left (+ strand)
b.add_patch(mp.Rectangle((tss - 1.0, 0.95), gx1 - tss + 1.0, 3.65, fc=tint(0.07), ec="none", zorder=0))
b.text(tss - 0.95, 4.68, "promoter + gene body", ha="left", va="bottom", fontsize=8.2, color=CRIMSON)
b.plot([0.9, 9.7], [4.0, 4.0], color=LINE, lw=1.0)             # genome axis
b.plot([tss, gx1], [4.0, 4.0], color=DARK, lw=1.6)
for ex0, ex1 in ((tss, tss + 0.45), (5.0, 5.35), (6.4, 6.8), (8.1, gx1)):
    b.add_patch(mp.Rectangle((ex0, 3.85), ex1 - ex0, 0.3, fc=DARK, ec="none", zorder=2))
b.annotate("", (tss + 0.55, 4.45), (tss, 4.2), arrowprops=dict(arrowstyle="-|>", color=DARK, lw=1.0,
                                                               connectionstyle="angle,angleA=90,angleB=0,rad=0",
                                                               mutation_scale=8))
b.text(gx1, 3.6, "gene G", ha="right", va="center", fontsize=8.6, color=TEXT, fontweight="semibold")
rows = [("cell 1", [2.95, 3.05, 3.2, 3.5, 3.62, 5.1, 6.5, 8.3, 1.4]),
        ("cell 2", [3.1, 3.55, 7.2, 9.3]),
        ("cell 3", [1.5, 9.4])]
for i, (nm, cuts) in enumerate(rows):
    y = 2.95 - i * 0.78
    b.plot([0.9, 9.7], [y, y], color=GREYF, lw=0.8)
    b.text(0.75, y, nm, ha="right", va="center", fontsize=8.6, color=TEXT2)
    inside = [c for c in cuts if tss - 1.0 <= c <= gx1]
    for c in cuts:
        b.plot([c, c], [y - 0.17, y + 0.17], color=CRIMSON if c in inside else "#9AA1AB", lw=1.6)
    b.text(9.95, y, f"{len(inside)}", ha="center", va="center", fontsize=10, color=CRIMSON, fontweight="bold")
b.text(9.95, 3.5, "sum", ha="center", va="center", fontsize=8.2, color=TEXT2)
b.text(0.2, 0.6, "Count each cell's reads over the gene's promoter and body → a gene-like score.\n"
       "Cells with a high score for marker genes can borrow RNA cell-type names.", ha="left", va="center",
       fontsize=8.6, color=TEXT, linespacing=1.35)
b.set_title("(b) Gene activity score", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=24)
b.text(0, 1.03, "ticks = Tn5 cut sites read in each cell", transform=b.transAxes, ha="left", va="bottom", fontsize=9.1,
       color=TEXT2)
header(axes[0], "Beyond scRNA-seq · single-cell chromatin", "snATAC-seq: open chromatin, one nucleus at a time",
       dy=44, x1=2.21)
fig.savefig(OUT + "snatac.png"); print("ok k34 snatac", (M == 0).mean().round(2))
