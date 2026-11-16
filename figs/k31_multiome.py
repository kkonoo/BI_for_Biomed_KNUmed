from gallery import *
rng = np.random.default_rng(7)

# snATAC-seq gives a very sparse cells × peaks matrix; multiome reads RNA and ATAC from the same nucleus,
# so the shared barcode links a peak to a gene cell by cell
fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.3), gridspec_kw={"wspace": 0.1})
DARK = "#5E6670"; GREYF = "#E3E6EB"
for ax in axes:
    ax.set_xlim(0, 10); ax.set_ylim(0, 6); ax.axis("off")

def cell(ax, x, y, v, s=0.36):
    ax.add_patch(mp.Rectangle((x, y), s, s, fc=tint(v) if v > 0 else GREYF, ec="white", lw=1.2))

# (a) snATAC-seq: cells × peaks, mostly 0
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
a.set_title("(a) snATAC-seq", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=24)
a.text(0, 1.03, "open chromatin, one nucleus at a time", transform=a.transAxes, ha="left", va="bottom",
       fontsize=9.1, color=TEXT2)

# (b) multiome: same barcode → RNA and ATAC per cell → peak–gene link
b = axes[1]; n = 8; s = 0.5; x0 = 1.9
g = np.array([0.95, 0.85, 0.7, 0.1, 0.05, 0.6, 0.0, 0.9])
p = np.array([1.0, 1.0, 0.55, 0.0, 0.0, 0.55, 0.0, 1.0])
for j in range(n):
    b.text(x0 + j * s + s / 2, 4.75, f"{j + 1}", ha="center", va="center", fontsize=8.4, color=TEXT2)
    cell(b, x0 + j * s + 0.04, 3.85, g[j], s=s - 0.08)
    cell(b, x0 + j * s + 0.04, 2.95, p[j], s=s - 0.08)
b.text(x0 + n * s / 2, 5.2, "cells (one barcode each)", ha="center", fontsize=8.8, color=TEXT2)
b.text(x0 - 0.15, 4.06, "RNA: gene G", ha="right", va="center", fontsize=8.8, color=TEXT)
b.text(x0 - 0.15, 3.16, "ATAC: peak P", ha="right", va="center", fontsize=8.8, color=TEXT)
b.annotate("", (x0 + n * s + 0.35, 3.16), (x0 + n * s + 0.35, 4.06),
           arrowprops=dict(arrowstyle="<->", color=CRIMSON, lw=1.2))
b.text(x0 + n * s + 0.5, 3.61, "move\ntogether", ha="left", va="center", fontsize=8.6, color=CRIMSON,
       fontweight="semibold", linespacing=1.2)
b.text(0.4, 1.75, "Same nucleus, same barcode → RNA and open chromatin\nare paired cell by cell", ha="left",
       va="center", fontsize=8.8, color=TEXT, linespacing=1.35)
b.text(0.4, 0.75, "→ peak P opens where gene G is on: a candidate enhancer\n   (a hypothesis to test, as in W3)",
       ha="left", va="center", fontsize=8.6, color=CRIMSON, linespacing=1.35)
b.set_title("(b) Multiome (RNA + ATAC)", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=24)
b.text(0, 1.03, "two readouts from the same cell", transform=b.transAxes, ha="left", va="bottom", fontsize=9.1,
       color=TEXT2)

header(axes[0], "Beyond scRNA-seq · single-cell chromatin", "Reading chromatin in single cells, alone or with RNA",
       dy=44, x1=2.1)
fig.savefig(OUT + "multiome.png"); print("ok k31 multiome")
