from gallery import *
rng = np.random.default_rng(14)

# Biological vs technical replicates: only independent individuals count toward n
fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.9), gridspec_kw={"wspace": 0.1})
DARK = "#5E6670"; GREYF = "#B8BEC7"; BODY = "#D5D9DF"

def mouse(ax, x, y, c=BODY):
    ax.add_patch(mp.Ellipse((x, y), 0.9, 0.5, fc=c, ec="none"))
    ax.add_patch(mp.Circle((x + 0.5, y + 0.12), 0.2, fc=c, ec="none"))
    ax.add_patch(mp.Circle((x + 0.5, y + 0.34), 0.09, fc=c, ec="none"))
    ax.plot([x - 0.45, x - 0.75, x - 0.85], [y, y + 0.05, y + 0.2], color=c, lw=1.6)
def tube(ax, x, y, c):
    ax.add_patch(mp.FancyBboxPatch((x - 0.12, y - 0.32), 0.24, 0.64, boxstyle="round,pad=0,rounding_size=0.1",
                                   fc="white", ec=TEXT2, lw=1))
    ax.add_patch(mp.FancyBboxPatch((x - 0.1, y - 0.3), 0.2, 0.3, boxstyle="round,pad=0,rounding_size=0.08",
                                   fc=c, ec="none"))
def lib(ax, x, y, c):
    ax.add_patch(mp.Rectangle((x - 0.25, y - 0.12), 0.5, 0.24, fc=c, ec="none"))
def arrow(ax, a, b):
    ax.annotate("", b, a, arrowprops=dict(arrowstyle="-|>", color=LINE, lw=1, mutation_scale=8, shrinkA=2, shrinkB=2))

def panel(ax, title, sub, bio):
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
    ax.set_title(title, loc="left", fontsize=12.5, fontweight="semibold", color=CRIMSON if bio else TEXT, pad=22)
    ax.text(0, 1.03, sub, transform=ax.transAxes, ha="left", va="bottom", fontsize=9.2, color=TEXT2)
    cols = [tint(0.75), tint(0.45), tint(0.25)] if bio else [GREYF] * 3
    ys = [8.2, 6.9, 5.6]
    for k, y in enumerate(ys):
        if bio or k == 1:
            mouse(ax, 1.3, y, BODY)
            ax.text(1.3, y - 0.55, f"mouse {k + 1}" if bio else "mouse 1", ha="center", va="center", fontsize=8, color=TEXT2)
        if bio or k == 1:
            tube(ax, 3.6, y, cols[k])
            arrow(ax, (2.1, y), (3.4, y))
        src = (3.8, y) if bio else (3.8, 6.9)
        lib(ax, 6.0, y, cols[k]); arrow(ax, src, (5.7, y))
    ax.text(3.6, 9.25, "RNA", ha="center", va="center", fontsize=8.5, color=TEXT2)
    ax.text(6.0, 9.25, "library →\nsequencing", ha="center", va="center", fontsize=8.5, color=TEXT2, linespacing=1.15)
    # measured expression of one gene
    ax.plot([0.3, 9.7], [2.6, 2.6], color=LINE, lw=1)
    vals = [3.0, 5.4, 7.6] if bio else [4.9, 5.2, 5.5]
    for v, c in zip(vals, cols):
        ax.plot([v], [2.6], "o", ms=10 if bio else 8, color=c, mec="white", mew=1.2)
    ax.text(0.3, 3.35, "expression of one gene", ha="left", va="center", fontsize=8.8, color=TEXT2)
    ax.text(5.0, 1.75, "spread = real variation between individuals" if bio else "almost identical: same RNA, read 3×",
            ha="center", va="center", fontsize=9, color=TEXT2)
    ax.text(5.0, 0.55, "n = 3" if bio else "n = 1, not 3", ha="center", va="center", fontsize=14, fontweight="bold",
            color=CRIMSON if bio else DARK)

panel(axes[0], "Biological replicates", "different animals, patients, or independent cultures", True)
panel(axes[1], "Technical replicates", "the same RNA sequenced more than once", False)
fig.text(0.5, -0.02, "Three flasks split from one culture on the same day are closer to technical than biological replicates",
         ha="center", va="center", fontsize=9, color=TEXT2)
header(axes[0], "Bulk RNA-seq · replicates", "Only biological replicates count toward n", dy=40, x1=2.1)
fig.savefig(OUT + "replicates.png"); print("ok h14 replicates")
