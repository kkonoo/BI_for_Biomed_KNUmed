from gallery import *

# Microwell loading (BD Rhapsody): cells settle by gravity into wells, then one barcoded bead drops into each well,
# cells are lysed in place and their mRNA is caught by the bead in the same well; beads are then pooled
W, H = 9.8, 3.55
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, (H - 0.75) / H]); ax.set_xlim(0, W); ax.set_ylim(0, H - 0.75); ax.axis("off")
ax.set_aspect("equal")
DARK = "#5E6670"; PLATE = "#D5D9DF"; BEAD = tint(0.45)
WW, WD, STEP, N = 0.36, 0.6, 0.45, 6          # well width, depth, pitch, wells per panel
YB = 1.3                                         # plate top


def plate(x0):
    """Cross-section of a microwell plate: a grey slab with N open wells; returns the well centres."""
    xs = [x0 + 0.12 + WW / 2 + i * STEP for i in range(N)]
    x1 = xs[-1] + WW / 2 + 0.12
    ax.add_patch(mp.Rectangle((x0, YB - WD - 0.18), x1 - x0, WD + 0.18, fc=PLATE, ec="none", zorder=1))
    for x in xs:
        ax.add_patch(mp.Rectangle((x - WW / 2, YB - WD), WW, WD + 0.01, fc="white", ec="none", zorder=2))
    return xs


def squiggle(x, y, c=CRIMSON, s=0.09):
    t = np.linspace(0, 1, 30)
    ax.plot(x + s * 1.6 * t, y + s * 0.35 * np.sin(t * 4 * np.pi), color=c, lw=1.0, zorder=5)


panels = [(0.1, "①  Cells settle by gravity", "loaded dilute → most wells stay empty"),
          (3.45, "②  One bead per well", "bead ≈ well size → only one fits"),
          (6.8, "③  Lyse and capture", "mRNA sticks to the bead in its own well")]
cells = [0, 1, 0, 0, 2, 0]
for k, (x0, t, s) in enumerate(panels):
    ax.text(x0, 2.57, t, ha="left", va="center", fontsize=11.2, fontweight="semibold", color=TEXT)
    ax.text(x0, 2.3, s, ha="left", va="center", fontsize=8.6, color=TEXT2)
    xs = plate(x0)
    for x, n in zip(xs, cells):
        for j in range(n):                       # cells sit at the bottom of the well
            cx, cy = x - 0.068 * (n - 1) + 0.136 * j, YB - WD + 0.085
            if k < 2:
                ax.add_patch(mp.Circle((cx, cy), 0.065, fc=DARK, ec="none", zorder=3))
            else:
                ax.add_patch(mp.Circle((cx, cy), 0.065, fc="none", ec=DARK, lw=0.8, ls=(0, (2, 1.5)), zorder=3))
        if k >= 1:                               # one bead on top
            ax.add_patch(mp.Circle((x, YB - 0.2), 0.16, fc=BEAD, ec="none", zorder=4))
        if k == 2 and n:
            for dy in (-0.05, 0.06):
                squiggle(x - 0.08, YB - 0.2 + dy)
    if k == 0:                                   # cells still falling
        for x, y in ((xs[1] + 0.05, 1.95), (xs[3] - 0.05, 1.75), (xs[5] + 0.06, 2.0)):
            ax.add_patch(mp.Circle((x, y), 0.075, fc=DARK, ec="none", zorder=3))
            ax.annotate("", (x, y - 0.25), (x, y - 0.1), arrowprops=dict(arrowstyle="-|>", color=TEXT2, lw=0.8,
                                                                         mutation_scale=7))
    if k == 1:
        ax.text(xs[4], 0.42, "2 cells, one bead → a doublet", ha="center", va="top", fontsize=8, color=CRIMSON,
                linespacing=1.2)
    if k < 2:
        ax.annotate("", (panels[k + 1][0] - 0.06, YB - 0.3), (xs[-1] + WW / 2 + 0.2, YB - 0.3),
                    arrowprops=dict(arrowstyle="-|>", color=LINE, lw=1.3, mutation_scale=11))
ax.text(6.8, 0.42, "then the beads are pulled out with a magnet and pooled;\nthe barcode on each bead marks its well (cell)",
        ha="left", va="top", fontsize=8, color=TEXT2, linespacing=1.25)
header(ax, "Single-cell RNA-seq · microwells", "Microwells: the cell goes in first, then the bead", dy=4)
fig.savefig(OUT + "microwell_loading.png"); print("ok k32 microwell")
