from gallery import *

# Cell–cell communication: what the algorithm measured (two mRNAs high) vs what the figure claims (a signal)
W, H = 9.8, 4.7
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, (H - 0.8) / H]); ax.set_xlim(0, W); ax.set_ylim(0, H - 0.8); ax.axis("off")
DARK = "#5E6670"; GREYF = "#B8BEC7"
def box(x, y, w, h, ec=LINE, fc="white"):
    ax.add_patch(mp.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.08", fc=fc, ec=ec, lw=0.9))
def cell(x, y, name, gene, c):
    ax.add_patch(mp.Circle((x, y), 0.42, fc=c, ec="none"))
    ax.text(x, y, name, ha="center", va="center", fontsize=8.6, color="white", fontweight="semibold")
    ax.text(x, y - 0.55, gene + " mRNA ↑", ha="center", va="center", fontsize=8, color=TEXT, fontstyle="italic")
# what the algorithm has
box(0.1, 1.85, 4.5, 1.75)
ax.text(0.3, 3.38, "What the algorithm has", ha="left", va="center", fontsize=11, fontweight="semibold", color=TEXT)
cell(1.05, 2.5, "macrophage", "TNF", DARK); cell(3.65, 2.5, "T cell", "TNFRSF1A", GREYF)
ax.plot([1.5, 3.2], [2.5, 2.5], color=TEXT2, lw=1, ls=(0, (3, 2)))
ax.text(2.35, 2.62, "both mRNAs high\nin the same\ndissociated sample", ha="center", va="bottom", fontsize=8.2,
        color=TEXT2, linespacing=1.15)
ax.annotate("", (5.15, 2.72), (4.7, 2.72), arrowprops=dict(arrowstyle="-|>", color=LINE, lw=1.1, mutation_scale=9))
# what the figure says
box(5.25, 1.85, 4.45, 1.75, ec=CRIMSON)
ax.text(5.45, 3.38, "What the figure says", ha="left", va="center", fontsize=11, fontweight="semibold", color=CRIMSON)
cell(6.15, 2.5, "macrophage", "TNF", DARK); cell(8.85, 2.5, "T cell", "TNFRSF1A", GREYF)
ax.annotate("", (8.4, 2.5), (6.6, 2.5), arrowprops=dict(arrowstyle="-|>", color=CRIMSON, lw=2.6, mutation_scale=16))
ax.text(7.5, 2.68, "TNF signalling", ha="center", va="bottom", fontsize=8.8, color=CRIMSON, fontweight="semibold")
ax.text(7.5, 2.3, "p < 0.001", ha="center", va="top", fontsize=8.8, color=CRIMSON, fontweight="semibold")
# three gaps
gaps = [("① Never adjacent?", "dissociation threw away\nall spatial information", "→ spatial data (§7.3)"),
        ("② Protein never made?", "mRNA is not protein: secretion,\ncleavage and trafficking are unseen", "→ proteomics (W6)"),
        ("③ Database not neutral", "well-studied pairs are\nover-represented", "→ rediscovers the literature")]
for i, (h, d, fix) in enumerate(gaps):
    x = 0.1 + i * 3.23
    ax.text(x, 1.45, h, ha="left", va="center", fontsize=10, fontweight="semibold", color=TEXT)
    ax.text(x, 0.95, d, ha="left", va="center", fontsize=8.5, color=TEXT2, linespacing=1.25)
    ax.text(x, 0.45, fix, ha="left", va="center", fontsize=8.6, color=CRIMSON, fontweight="semibold")
ax.text(0.1, 0.0, "Use these tools to choose what to test (\"we will test whether…\"), not as a result (\"we showed…\").",
        ha="left", va="center", fontsize=8.8, color=TEXT)
header(ax, "Single-cell RNA-seq · cell–cell communication", "Co-expression is not interaction", dy=0)
fig.savefig(OUT + "ccc_fallacy.png"); print("ok k17 ccc")
