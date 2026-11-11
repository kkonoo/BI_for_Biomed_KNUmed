from gallery import *
rng = np.random.default_rng(13)

# Stranded vs unstranded: reads where two genes on opposite strands overlap
fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.4), gridspec_kw={"wspace": 0.08})
DARK = "#5E6670"; GREYF = "#B8BEC7"
P0, P1, Q0, Q1 = 0.6, 6.4, 3.6, 9.4                      # gene P on +, gene Q on −; overlap 3.6–6.4
xs = [3.8, 5.3, 4.3, 5.6, 3.9, 5.0, 4.6]
strand = [1, -1, 1, 1, -1, 1, -1]
lanes = [0, 0, 1, 1, 2, 2, 3]

def draw(ax, stranded, title, sub, note, ncol):
    ax.set_xlim(0, 10); ax.set_ylim(-0.6, 5.2); ax.axis("off")
    ax.add_patch(mp.Rectangle((Q0, -0.25), P1 - Q0, 4.5, fc=FAINT, ec="none", zorder=0))
    ax.text((Q0 + P1) / 2, 4.42, "overlap", ha="center", va="center", fontsize=8.5, color=TEXT2)
    y = 0.5
    ax.plot([0.2, 9.8], [y, y], color=LINE, lw=1)
    ax.annotate("", (P1, y + 0.22), (P0, y + 0.22), arrowprops=dict(arrowstyle="-|>", color=CRIMSON, lw=4,
                                                                    mutation_scale=11, shrinkA=0, shrinkB=0))
    ax.annotate("", (Q0, y - 0.22), (Q1, y - 0.22), arrowprops=dict(arrowstyle="-|>", color=DARK, lw=4,
                                                                    mutation_scale=11, shrinkA=0, shrinkB=0))
    ax.text(P0, y + 0.55, "gene P  (+ strand)", ha="left", va="center", fontsize=8.8, color=CRIMSON, fontweight="semibold")
    ax.text(Q1, y - 0.6, "gene Q  (− strand)", ha="right", va="center", fontsize=8.8, color=DARK, fontweight="semibold")
    for x, s, l in zip(xs, strand, lanes):
        yy = 1.45 + l * 0.6
        if stranded:
            col = CRIMSON if s > 0 else DARK
            a, b = ((x, yy), (x + 0.65, yy)) if s > 0 else ((x + 0.65, yy), (x, yy))
            ax.annotate("", b, a, arrowprops=dict(arrowstyle="-|>", color=col, lw=1.6, mutation_scale=8))
        else:
            ax.plot([x, x + 0.65], [yy, yy], color=GREYF, lw=2.6, solid_capstyle="round")
    ax.set_title(title, loc="left", fontsize=12.5, fontweight="semibold", color=ncol, pad=22)
    ax.text(0, 1.03, sub, transform=ax.transAxes, ha="left", va="bottom", fontsize=9.2, color=TEXT2)
    ax.text(5.0, -0.55, note, ha="center", va="center", fontsize=9.4, fontweight="semibold",
            color=CRIMSON if stranded else TEXT2)

draw(axes[0], False, "Unstranded", "reads do not record their strand",
     "P or Q?  → ambiguous: dropped or guessed", TEXT)
draw(axes[1], True, "Stranded", "each read keeps the strand it came from",
     "→ reads go to the right gene", CRIMSON)
header(axes[0], "Bulk RNA-seq · strandedness", "Where genes overlap, only a stranded library can tell them apart",
       dy=40, x1=2.08)
fig.savefig(OUT + "stranded.png"); print("ok h13 stranded")
