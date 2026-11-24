from gallery import *

# Multiome reads RNA and ATAC from the same nucleus: (a) both libraries carry the same cell barcode,
# (b) so a peak can be linked to a gene cell by cell (snATAC-seq alone is drawn in k34_snatac.py)
fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.3), gridspec_kw={"wspace": 0.1, "width_ratios": [1.12, 1]})
DARK = "#5E6670"; GREYF = "#E3E6EB"
for ax in axes:
    ax.set_xlim(0, 10); ax.set_ylim(0, 6); ax.axis("off")

def cell(ax, x, y, v, s=0.36):
    ax.add_patch(mp.Rectangle((x, y), s, s, fc=tint(v) if v > 0 else GREYF, ec="white", lw=1.2))

# (a) one nucleus, two libraries that carry the same cell barcode
a = axes[0]
a.add_patch(mp.Ellipse((1.25, 3.2), 2.0, 1.75, fc=FAINT, ec=LINE, lw=1.0))
t = np.linspace(0, 1, 60)
a.plot(0.75 + 1.0 * t, 3.45 + 0.18 * np.sin(t * 9 * np.pi), color=DARK, lw=1.4)          # chromatin
for dx, dy in ((0.0, 0.0), (0.35, -0.3)):                                                  # mRNA
    a.plot(0.85 + dx + 0.45 * t, 2.75 + dy + 0.06 * np.sin(t * 5 * np.pi), color=CRIMSON, lw=1.2)
a.text(1.25, 2.15, "one nucleus", ha="center", va="top", fontsize=8.6, color=TEXT2)
libs = [(4.35, "RNA library", [("barcode", DARK, 1.6), ("UMI", CRIMSON, 0.9), ("cDNA → gene", tint(0.3), 2.6)]),
        (2.25, "ATAC library", [("barcode", DARK, 1.6), ("Tn5 fragment → peak", "#C9CED6", 3.5)])]
for y, nm, parts in libs:
    a.annotate("", (3.05, y), (2.25, 3.2 + (y - 3.3) * 0.45), arrowprops=dict(arrowstyle="-|>", color=LINE, lw=1.2,
                                                                         mutation_scale=10))
    a.text(3.2, y + (0.45 if y > 3 else -0.45), nm, ha="left", va="center", fontsize=9.2, color=TEXT, fontweight="semibold")
    x = 3.2
    for lab, c, w in parts:
        a.add_patch(mp.Rectangle((x, y - 0.2), w, 0.4, fc=c, ec="white", lw=1.2))
        a.text(x + w / 2, y, lab, ha="center", va="center", fontsize=8.1, fontweight="semibold",
               color="white" if c in (DARK, CRIMSON) else TEXT)
        x += w
a.annotate("", (4.0, 2.47), (4.0, 4.13), arrowprops=dict(arrowstyle="<|-|>", color=DARK, lw=0.9, mutation_scale=8))
a.text(4.15, 3.3, "from the same bead:\nbarcodes matched one-to-one", ha="left", va="center", fontsize=8.6, linespacing=1.2,
       color=TEXT, fontweight="semibold")
a.text(0.4, 0.75, "Both libraries are sequenced separately, then joined by barcode →\n"
       "a cells × genes and a cells × peaks matrix with the same rows (cells)", ha="left", va="center",
       fontsize=8.6, color=TEXT, linespacing=1.35)
a.set_title("(a) One nucleus, two readouts", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=24)
a.text(0, 1.03, "both readouts are tagged by the same bead", transform=a.transAxes, ha="left", va="bottom",
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
b.text(0.4, 1.75, "Same nucleus, matched barcodes → RNA and open chromatin\nare paired cell by cell", ha="left",
       va="center", fontsize=8.8, color=TEXT, linespacing=1.35)
b.text(0.4, 0.75, "→ peak P opens where gene G is on: a candidate enhancer\n   (a hypothesis to test, as in W3)",
       ha="left", va="center", fontsize=8.6, color=CRIMSON, linespacing=1.35)
b.set_title("(b) Pairing RNA and ATAC per cell", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=24)
b.text(0, 1.03, "a gene and a peak that move together", transform=b.transAxes, ha="left", va="bottom", fontsize=9.1,
       color=TEXT2)

header(axes[0], "Beyond scRNA-seq · multiome", "Multiome: RNA and ATAC from the same nucleus",
       dy=44, x1=2.1)
fig.savefig(OUT + "multiome.png"); print("ok k31 multiome")
