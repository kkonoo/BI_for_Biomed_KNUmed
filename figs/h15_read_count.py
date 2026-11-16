from gallery import *

# From reads to counts, one figure per decision:
# spliced alignment, pseudo-alignment, transcript → gene, multi-mapping reads
DARK = "#5E6670"; GREYF = "#B8BEC7"; EX = "#5E6670"
EYEBROW = "Bulk RNA-seq · from reads to counts"

def new(title):
    fig = plt.figure(figsize=(6.6, 3.5))
    ax = fig.add_axes([0.02, 0.02, 0.96, 0.74]); ax.set_xlim(0, 10); ax.set_ylim(0.2, 5.4); ax.axis("off")
    header(ax, EYEBROW, title)
    return fig, ax

def exon(ax, x0, x1, y, c=EX, h=0.42):
    ax.add_patch(mp.Rectangle((x0, y - h / 2), x1 - x0, h, fc=c, ec="none", zorder=3))

# 1 — spliced alignment
fig, a = new("A read from mature mRNA jumps over the intron")
a.plot([0.3, 9.7], [2.0, 2.0], color=LINE, lw=1.4)
exon(a, 1.0, 3.6, 2.0); exon(a, 6.6, 9.0, 2.0)
a.text(2.3, 1.35, "exon 1", ha="center", fontsize=9.5, color=TEXT2); a.text(7.8, 1.35, "exon 2", ha="center", fontsize=9.5, color=TEXT2)
a.text(5.1, 1.35, "intron (in the genome)", ha="center", fontsize=9.5, color=TEXT2)
y = 3.2
a.plot([2.6, 3.6], [y, y], color=CRIMSON, lw=4, solid_capstyle="butt")
a.plot([6.6, 7.6], [y, y], color=CRIMSON, lw=4, solid_capstyle="butt")
a.plot([3.6, 6.6], [y, y], color=CRIMSON, lw=1, ls=(0, (3, 2)))
a.text(5.1, y + 0.35, "gap = intron", ha="center", va="bottom", fontsize=9.5, color=CRIMSON, fontweight="semibold")
a.text(5.1, 4.55, "read from mature mRNA (exon 1 + exon 2)", ha="center", fontsize=9.5, color=TEXT)
a.text(0.3, 0.5, "STAR, HISAT2 · keeps the position → splicing, fusions, new transcripts, browser",
       ha="left", fontsize=9, color=TEXT2)
fig.savefig(OUT + "spliced_alignment.png"); plt.close(fig)

# 2 — pseudo-alignment
fig, b = new("Pseudo-alignment asks which transcripts a read fits")
tx = [("T1", [(0.9, 2.6), (3.2, 4.8), (5.6, 7.2)], True), ("T2", [(0.9, 2.6), (3.2, 4.8), (7.8, 9.2)], True),
      ("T3", [(0.9, 2.6), (5.6, 7.2), (7.8, 9.2)], False)]
for k, (name, exs, ok) in enumerate(tx):
    yy = 3.3 - k * 0.9
    b.plot([exs[0][0], exs[-1][1]], [yy, yy], color=LINE, lw=1)
    for x0, x1 in exs: exon(b, x0, x1, yy, c=EX if ok else GREYF, h=0.34)
    b.text(0.55, yy, name, ha="right", va="center", fontsize=10, fontweight="semibold", color=TEXT)
    b.text(9.55, yy, "✓" if ok else "✕", ha="left", va="center", fontsize=11, color=CRIMSON if ok else DARK, fontweight="bold")
for i in range(4):                                            # the read as k-mers spanning exon 1 → exon 2
    b.add_patch(mp.Rectangle((2.0 + i * 0.55, 4.45), 0.5, 0.28, fc=tint(0.55), ec="white", lw=1))
b.text(4.4, 4.59, "read → k-mers (exon 1 | exon 2)", ha="left", va="center", fontsize=9.3, color=TEXT)
b.text(0.3, 0.5, "Salmon, kallisto · compatible with {T1, T2} · no base-by-base alignment",
       ha="left", fontsize=9, color=TEXT2)
fig.savefig(OUT + "pseudo_alignment.png"); plt.close(fig)

# 3 — transcript → gene
fig, c = new("Estimate each isoform, then sum to the gene")
for k, (name, val, col) in enumerate((("isoform 1", 60, tint(0.6)), ("isoform 2", 40, tint(0.3)))):
    yy = 4.3 - k * 1.0
    c.text(0.3, yy, name, ha="left", va="center", fontsize=10, color=TEXT)
    c.add_patch(mp.Rectangle((2.2, yy - 0.25), val / 100 * 5.5, 0.5, fc=col, ec="none"))
    c.text(2.3 + val / 100 * 5.5, yy, f"{val}", ha="left", va="center", fontsize=10, color=TEXT2)
c.plot([2.2, 7.7], [2.55, 2.55], color=LINE, lw=1)
c.text(0.3, 1.8, "gene", ha="left", va="center", fontsize=10.5, fontweight="bold", color=TEXT)
c.add_patch(mp.Rectangle((2.2, 1.55), 3.3, 0.5, fc=tint(0.6), ec="none"))
c.add_patch(mp.Rectangle((5.5, 1.55), 2.2, 0.5, fc=tint(0.3), ec="none"))
c.text(7.85, 1.8, "100", ha="left", va="center", fontsize=10.5, fontweight="bold", color=CRIMSON)
c.text(0.3, 0.5, "tximport · the isoform split is noisy; the gene total is robust", ha="left", fontsize=9,
       color=TEXT2)
fig.savefig(OUT + "tx_to_gene.png"); plt.close(fig)

# 4 — multi-mapping
fig, d = new("One read, two equally good places")
for k, (name, yy) in enumerate((("gene A", 3.9), ("paralogue A′", 2.4))):
    d.plot([0.6, 9.4], [yy, yy], color=LINE, lw=1)
    exon(d, 3.0, 7.0, yy, c=EX, h=0.34)
    d.text(0.6, yy + 0.4, name, ha="left", va="center", fontsize=9.7, color=TEXT)
    d.plot([4.4, 5.6], [yy + 0.45, yy + 0.45], color=CRIMSON, lw=3.5, solid_capstyle="butt")
d.annotate("", (5.0, 2.9), (5.0, 3.4), arrowprops=dict(arrowstyle="<->", color=CRIMSON, lw=1))
d.text(5.25, 3.15, "same read, two equal hits", ha="left", va="center", fontsize=9.3, color=CRIMSON)
d.text(0.3, 1.25, "drop it → gene families are under-counted", ha="left", fontsize=9.3, color=TEXT)
d.text(0.3, 0.6, "split it → better, but adds its own uncertainty", ha="left", fontsize=9.3, color=TEXT)
fig.savefig(OUT + "multimapping.png"); plt.close(fig)
print("ok h15 read to count (4 figures)")
