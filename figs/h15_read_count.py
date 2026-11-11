from gallery import *

# From reads to counts: spliced alignment, pseudo-alignment, transcript → gene, multi-mapping reads
fig, axs = plt.subplots(2, 2, figsize=(9.8, 7.4), gridspec_kw={"wspace": 0.12, "hspace": 0.5})
DARK = "#5E6670"; GREYF = "#B8BEC7"; EX = "#5E6670"
def panel(ax, t, sub):
    ax.set_xlim(0, 10); ax.set_ylim(0, 6); ax.axis("off")
    ax.set_title(t, loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
    ax.text(0, 1.03, sub, transform=ax.transAxes, ha="left", va="bottom", fontsize=9.1, color=TEXT2)
def exon(ax, x0, x1, y, c=EX, h=0.42):
    ax.add_patch(mp.Rectangle((x0, y - h / 2), x1 - x0, h, fc=c, ec="none", zorder=3))

# (a) spliced alignment
a = axs[0, 0]; panel(a, "(a) Spliced alignment", "the read jumps over the intron (STAR, HISAT2)")
a.plot([0.3, 9.7], [2.0, 2.0], color=LINE, lw=1.4)
exon(a, 1.0, 3.6, 2.0); exon(a, 6.6, 9.0, 2.0)
a.text(2.3, 1.35, "exon 1", ha="center", fontsize=9, color=TEXT2); a.text(7.8, 1.35, "exon 2", ha="center", fontsize=9, color=TEXT2)
a.text(5.1, 1.35, "intron (in the genome)", ha="center", fontsize=9, color=TEXT2)
y = 3.2
a.plot([2.6, 3.6], [y, y], color=CRIMSON, lw=4, solid_capstyle="butt")
a.plot([6.6, 7.6], [y, y], color=CRIMSON, lw=4, solid_capstyle="butt")
a.plot([3.6, 6.6], [y, y], color=CRIMSON, lw=1, ls=(0, (3, 2)))
a.text(5.1, y + 0.35, "gap = intron", ha="center", va="bottom", fontsize=9, color=CRIMSON, fontweight="semibold")
a.text(5.1, 4.55, "read from mature mRNA (exon 1 + exon 2)", ha="center", fontsize=9, color=TEXT)
a.text(0.3, 0.45, "keeps the position → splicing, fusions, new transcripts, browser", ha="left", fontsize=8.8, color=TEXT2)

# (b) pseudo-alignment
b = axs[0, 1]; panel(b, "(b) Pseudo-alignment", "which transcripts could the read come from? (Salmon, kallisto)")
tx = [("T1", [(0.9, 2.6), (3.2, 4.8), (5.6, 7.2)], True), ("T2", [(0.9, 2.6), (3.2, 4.8), (7.8, 9.2)], True),
      ("T3", [(0.9, 2.6), (5.6, 7.2), (7.8, 9.2)], False)]
for k, (name, exs, ok) in enumerate(tx):
    yy = 3.3 - k * 0.9
    b.plot([exs[0][0], exs[-1][1]], [yy, yy], color=LINE, lw=1)
    for x0, x1 in exs: exon(b, x0, x1, yy, c=EX if ok else GREYF, h=0.34)
    b.text(0.55, yy, name, ha="right", va="center", fontsize=9.5, fontweight="semibold", color=TEXT)
    b.text(9.55, yy, "✓" if ok else "✕", ha="left", va="center", fontsize=11, color=CRIMSON if ok else DARK, fontweight="bold")
for i in range(4):                                            # the read as k-mers spanning exon 1 → exon 2
    b.add_patch(mp.Rectangle((2.0 + i * 0.55, 4.45), 0.5, 0.28, fc=tint(0.55), ec="white", lw=1))
b.text(2.0, 5.0, "read → k-mers (exon 1 | exon 2)", ha="left", fontsize=8.8, color=TEXT)
b.text(0.3, 0.45, "compatible with {T1, T2} · no base-by-base alignment", ha="left", fontsize=8.8, color=TEXT2)

# (c) transcript → gene
c = axs[1, 0]; panel(c, "(c) Transcript level → gene level", "estimate per isoform, then sum (tximport)")
for k, (name, val, col) in enumerate((("isoform 1", 60, tint(0.6)), ("isoform 2", 40, tint(0.3)))):
    yy = 4.3 - k * 1.0
    c.text(0.3, yy, name, ha="left", va="center", fontsize=9.5, color=TEXT)
    c.add_patch(mp.Rectangle((2.2, yy - 0.25), val / 100 * 5.5, 0.5, fc=col, ec="none"))
    c.text(2.3 + val / 100 * 5.5, yy, f"{val}", ha="left", va="center", fontsize=9.5, color=TEXT2)
c.plot([2.2, 7.7], [2.55, 2.55], color=LINE, lw=1)
c.text(0.3, 1.8, "gene", ha="left", va="center", fontsize=10, fontweight="bold", color=TEXT)
c.add_patch(mp.Rectangle((2.2, 1.55), 3.3, 0.5, fc=tint(0.6), ec="none"))
c.add_patch(mp.Rectangle((5.5, 1.55), 2.2, 0.5, fc=tint(0.3), ec="none"))
c.text(7.85, 1.8, "100", ha="left", va="center", fontsize=10, fontweight="bold", color=CRIMSON)
c.text(0.3, 0.45, "isoform split is noisy; the gene total is robust", ha="left", fontsize=8.8, color=TEXT2)

# (d) multi-mapping
d = axs[1, 1]; panel(d, "(d) Multi-mapping reads", "fits two places equally well (paralogues, repeats)")
for k, (name, yy) in enumerate((("gene A", 3.6), ("paralogue A′", 2.1))):
    d.plot([0.6, 9.4], [yy, yy], color=LINE, lw=1)
    exon(d, 3.0, 7.0, yy, c=EX, h=0.34)
    d.text(0.6, yy + 0.4, name, ha="left", va="center", fontsize=9.2, color=TEXT)
    d.plot([4.4, 5.6], [yy + 0.45, yy + 0.45], color=CRIMSON, lw=3.5, solid_capstyle="butt")
d.annotate("", (5.0, 2.6), (5.0, 3.1), arrowprops=dict(arrowstyle="<->", color=CRIMSON, lw=1))
d.text(5.25, 2.85, "same read, two equal hits", ha="left", va="center", fontsize=8.8, color=CRIMSON)
d.text(0.3, 1.05, "drop it → gene families are under-counted", ha="left", fontsize=8.8, color=TEXT)
d.text(0.3, 0.45, "split it → better, but adds its own uncertainty", ha="left", fontsize=8.8, color=TEXT)

header(axs[0, 0], "Bulk RNA-seq · from reads to counts", "Every count matrix hides these four decisions",
       dy=40, x1=2.12)
fig.savefig(OUT + "read_to_count.png"); print("ok h15 read to count")
