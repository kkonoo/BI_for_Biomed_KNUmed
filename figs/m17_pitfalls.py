from gallery import *

# Where proteomics errors enter, and the sign each one leaves in the data
W, H = 9.8, 4.5
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, (H - 0.8) / H]); ax.set_xlim(0, W); ax.set_ylim(0, H - 0.8); ax.axis("off")
DARK = "#5E6670"; GREYF = "#B8BEC7"
steps = [("Collection", "tube type, time\nto freezing", "red-cell and platelet\nproteins differ by group", True),
         ("Storage", "freeze–thaw\ncycles", "degradation differs\nby biobank or site", True),
         ("Preparation", "keratin from skin\nand hair; uneven\ndigestion", "keratins among\nthe top hits", True),
         ("Plasma\nhandling", "no depletion or\nfractionation", "mostly the most\nabundant proteins", True),
         ("Acquisition", "DDA sampling;\nTMT sets", "missingness or\nsignal tracks the batch", False),
         ("Analysis", "imputation;\nprotein groups", "big fold changes at\nthe detection limit", False)]
N = len(steps); CW, GAP = 1.42, 0.13
X0 = (W - (N * CW + (N - 1) * GAP)) / 2
for i, (nm, err, sign, bench) in enumerate(steps):
    x = X0 + i * (CW + GAP)
    col = CRIMSON if bench else DARK
    ax.add_patch(mp.FancyBboxPatch((x, 2.55), CW, 0.62, boxstyle="round,pad=0,rounding_size=0.08",
                                   fc=col, ec="none"))
    ax.text(x + CW / 2, 2.86, nm, ha="center", va="center", fontsize=9.6, fontweight="bold", color="white",
            linespacing=1.05)
    if i < N - 1:
        ax.annotate("", (x + CW + GAP, 2.86), (x + CW, 2.86),
                    arrowprops=dict(arrowstyle="-|>", color=GREYF, lw=1, mutation_scale=7, shrinkA=0, shrinkB=0))
    ax.text(x + 0.06, 2.3, "what goes wrong", ha="left", va="center", fontsize=7.6, color=TEXT2)
    ax.text(x + 0.06, 2.12, err, ha="left", va="top", fontsize=8.6, color=TEXT, linespacing=1.15)
    ax.plot([x + 0.06, x + CW - 0.06], [1.3, 1.3], color=LINE, lw=0.8)
    ax.text(x + 0.06, 1.12, "sign in the data", ha="left", va="center", fontsize=7.6, color=TEXT2)
    ax.text(x + 0.06, 0.94, sign, ha="left", va="top", fontsize=8.6, color=col, linespacing=1.15,
            fontweight="semibold")
xb = X0 + 4 * (CW + GAP) - GAP
ax.plot([X0, xb], [3.27, 3.27], color=CRIMSON, lw=1.4)
ax.text(X0, 3.36, "at the bench, before any software runs", ha="left", va="bottom", fontsize=9, color=CRIMSON,
        fontweight="semibold")
ax.plot([xb + GAP, W - X0], [3.27, 3.27], color=DARK, lw=1.4)
ax.text(xb + GAP, 3.36, "in the instrument and code", ha="left", va="bottom", fontsize=9, color=DARK,
        fontweight="semibold")
ax.text(X0, 0.12, "If one group was collected or stored differently, every protein that reacts to handling becomes "
        "a 'difference' between groups.", ha="left", va="center", fontsize=8.8, color=TEXT2)
header(ax, "Proteomics · pitfalls", "Most proteomics errors start at the bench, not in the software", x=X0 / W,
       dy=0, x1=1 - X0 / W)
fig.savefig(OUT + "pitfalls.png"); print("ok m17 pitfalls")
