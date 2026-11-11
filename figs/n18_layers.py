from gallery import *

# Beyond bacterial DNA: each layer answers a different question, with butyrate as the running example
W, H = 9.8, 4.5
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, (H - 0.8) / H]); ax.set_xlim(0, W); ax.set_ylim(0, H - 0.8); ax.axis("off")
DARK = "#5E6670"; GREYF = "#B8BEC7"
layers = [("DNA", "metagenome", "who is there, and\nwhich genes they carry", "butyrate-making\ngenes are present",
           "a gene present may\nnever be used", tint(0.35)),
          ("RNA", "metatranscriptome", "which genes are\nswitched on", "those genes are\nbeing expressed",
           "RNA decays fast;\nhard to sample well", tint(0.6)),
          ("small molecules", "metabolome", "what is actually made\nand reaches the host", "butyrate measured\nin stool or blood",
           "not every molecule\ncomes from microbes", CRIMSON)]
N = len(layers); CW, GAP = 2.45, 0.42
X0 = 0.2
for i, (nm, assay, q, ex, blind, col) in enumerate(layers):
    x = X0 + i * (CW + GAP)
    ax.add_patch(mp.FancyBboxPatch((x, 2.8), CW, 0.62, boxstyle="round,pad=0,rounding_size=0.08", fc=col, ec="none"))
    ax.text(x + 0.12, 3.22, nm, ha="left", va="center", fontsize=10.5, fontweight="bold",
            color="white" if col == CRIMSON else TEXT)
    ax.text(x + 0.12, 2.97, assay, ha="left", va="center", fontsize=8.6, color="white" if col == CRIMSON else TEXT)
    if i < N - 1:
        ax.annotate("", (x + CW + GAP - 0.04, 3.11), (x + CW + 0.04, 3.11),
                    arrowprops=dict(arrowstyle="-|>", color=GREYF, lw=1.1, mutation_scale=9, shrinkA=0, shrinkB=0))
    for k, (lab, txt, c, wt) in enumerate((("answers", q, TEXT, "normal"), ("butyrate example", ex, DARK, "semibold"),
                                           ("blind spot", blind, CRIMSON, "normal"))):
        y = 2.55 - k * 0.8
        ax.text(x + 0.05, y, lab, ha="left", va="center", fontsize=7.8, color=TEXT2)
        ax.text(x + 0.05, y - 0.16, txt, ha="left", va="top", fontsize=8.8, color=c, fontweight=wt, linespacing=1.15)
# host at the end
xh = X0 + N * (CW + GAP) - GAP + 0.25
ax.annotate("", (xh + 0.02, 3.11), (xh - 0.25, 3.11),
            arrowprops=dict(arrowstyle="-|>", color=GREYF, lw=1.1, mutation_scale=9, shrinkA=0, shrinkB=0))
ax.text(xh + 0.06, 3.11, "host", ha="left", va="center", fontsize=10.5, fontweight="bold", color=TEXT)
ax.text(X0, 0.08, "Not only bacteria: viruses (mostly phage) and fungi live there too, and most of them are missing "
        "from reference databases.", ha="left", va="center", fontsize=8.6, color=TEXT2)
header(ax, "Metagenomics · beyond DNA", "Present, active, and reaching the host are three different claims",
       x=X0 / W, dy=0, x1=1 - X0 / W)
fig.savefig(OUT + "omics_layers.png"); print("ok n18")
