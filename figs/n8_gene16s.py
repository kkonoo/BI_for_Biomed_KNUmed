from gallery import *

# The 16S rRNA gene: (a) variable regions sit between conserved stretches where primers bind (E. coli numbering);
# (b) equal numbers of cells do not give equal numbers of reads: copy number and primer fit both distort
fig = plt.figure(figsize=(9.8, 5.2))
gs = fig.add_gridspec(2, 2, height_ratios=[0.8, 1.25], width_ratios=[1.35, 1], hspace=0.62, wspace=0.3,
                      left=0.06, right=0.98, top=0.82, bottom=0.08)
DARK = "#5E6670"; GREYF = "#B8BEC7"

def ptitle(ax, t, sub):
    ax.set_title(t, loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
    ax.text(0, 1.03, sub, transform=ax.transAxes, fontsize=9, color=TEXT2)

a = fig.add_subplot(gs[0, :]); a.set_xlim(-40, 1560); a.set_ylim(-3.0, 2.1); a.axis("off")
L = 1542
V = [(69, 99), (137, 242), (433, 497), (576, 682), (822, 879), (986, 1043), (1117, 1173), (1243, 1294), (1435, 1465)]
a.add_patch(mp.Rectangle((0, 0), L, 0.8, fc="#E3E6EB", ec="none"))
for i, (s, e) in enumerate(V):
    a.add_patch(mp.Rectangle((s, 0), e - s, 0.8, fc=DARK, ec="none"))
    a.text((s + e) / 2, 1.0, f"V{i + 1}", ha="center", va="bottom", fontsize=8.6, color=TEXT)
a.text(L + 15, 0.4, "~1,540 bp", ha="left", va="center", fontsize=8.3, color=TEXT2)
pairs = [("V1–V2", 8, 355, "27F – 338R"), ("V3–V4", 341, 805, "341F – 805R"), ("V4", 515, 806, "515F – 806R")]
for k, (nm, s, e, pr) in enumerate(pairs):
    y = -0.75 - k * 0.85
    c = CRIMSON if nm == "V4" else tint(0.55)
    a.plot([s, e], [y, y], color=c, lw=1.2)
    a.annotate("", (s + 60, y), (s, y), arrowprops=dict(arrowstyle="-|>", color=c, lw=2, mutation_scale=9, shrinkA=0,
                                                         shrinkB=0))
    a.annotate("", (e - 60, y), (e, y), arrowprops=dict(arrowstyle="-|>", color=c, lw=2, mutation_scale=9, shrinkA=0,
                                                         shrinkB=0))
    a.text(e + 25, y, f"{nm}  ({pr})", ha="left", va="center", fontsize=8.4, color=TEXT)
a.text(1560, -1.6, "arrows = primers, which must bind\nconserved stretches (light grey)", ha="right", va="center",
       fontsize=8.4, color=TEXT2, linespacing=1.2)
ptitle(a, "(a) The gene and three common primer pairs", "dark = variable regions that tell taxa apart "
       "(E. coli numbering)")

# (b) three taxa, equal cells
b = fig.add_subplot(gs[1, 0])
copies = np.array([1, 4, 10]); eff = np.array([1.0, 1.0, 0.3])
cells = np.ones(3) / 3
genes = cells * copies; genes /= genes.sum()
reads = cells * copies * eff; reads /= reads.sum()
cols = [CRIMSON, DARK, GREYF]
names = ["taxon 1", "taxon 2", "taxon 3"]
stages = [("cells in\nthe sample", cells), ("16S genes\nin the tube", genes), ("reads after\nPCR", reads)]
for x, (nm, v) in enumerate(stages):
    bot = 0
    for k in range(3):
        b.bar(x, v[k], bottom=bot, color=cols[k], width=0.58, edgecolor="white", lw=1)
        if v[k] >= 0.09:
            b.text(x, bot + v[k] / 2, f"{v[k]:.0%}", ha="center", va="center", fontsize=8.6, fontweight="semibold",
                   color="white" if cols[k] != GREYF else TEXT)
        else:
            b.text(x + 0.33, bot + v[k] / 2, f"{v[k]:.0%}", ha="left", va="center", fontsize=8.6,
                   fontweight="semibold", color=cols[k])
        bot += v[k]
b.set_xticks(range(3)); b.set_xticklabels([s[0] for s in stages], fontsize=8.8, color=TEXT)
b.set_yticks([]); b.tick_params(axis="x", length=0)
for sp in ("left", "top", "right"): b.spines[sp].set_visible(False)
b.spines["bottom"].set_color(LINE)
b.set_xlim(-0.5, 2.5)
ptitle(b, "(b) Equal cells, unequal reads", "three taxa, one third of the cells each (computed)")

c = fig.add_subplot(gs[1, 1]); c.axis("off"); c.set_xlim(0, 1); c.set_ylim(0, 1)
rows = [("taxon 1", "1 copy per genome", "perfect primer match", CRIMSON),
        ("taxon 2", "4 copies", "perfect match", DARK),
        ("taxon 3", "10 copies", "one mismatch → 30% efficiency", GREYF)]
for k, (nm, cp, pm, col) in enumerate(rows):
    y = 0.9 - k * 0.27
    c.add_patch(mp.Rectangle((0, y - 0.05), 0.05, 0.1, fc=col, ec="none", transform=c.transAxes))
    c.text(0.09, y + 0.02, nm, ha="left", va="center", fontsize=9.5, fontweight="semibold", color=TEXT)
    c.text(0.09, y - 0.09, f"{cp}; {pm}", ha="left", va="center", fontsize=8.5, color=TEXT2)
c.text(0, 0.08, "16S copy number ranges from 1 to ~15 per\ngenome → read share ≠ cell share",
       ha="left", va="center", fontsize=8.8, color=CRIMSON, fontweight="semibold", linespacing=1.2)
header(a, "Metagenomics · 16S rRNA gene", "Conserved stretches let you amplify; variable regions let you classify",
       dy=40, x1=1)
fig.savefig(OUT + "gene_16s.png"); print("ok n8", np.round(genes, 3), np.round(reads, 3))
