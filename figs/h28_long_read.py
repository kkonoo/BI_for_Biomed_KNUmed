from gallery import *

# Why long reads resolve isoforms and fusions: short reads see pieces, a long read sees the whole molecule
fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.5), gridspec_kw={"wspace": 0.08})
DARK = "#5E6670"; GREYF = "#C3C8D0"
for ax in axes:
    ax.set_xlim(0, 10); ax.set_ylim(0, 6.2); ax.axis("off")

def box(ax, x0, x1, y, fc, h=0.36, label=None, tc="white"):
    ax.add_patch(mp.Rectangle((x0, y - h / 2), x1 - x0, h, fc=fc, ec="none", zorder=3))
    if label:
        ax.text((x0 + x1) / 2, y, label, ha="center", va="center", fontsize=7.6, color=tc, zorder=4,
                fontweight="semibold")

def read(ax, segs, y, col, lw=3.2):
    """A spliced read: thick bars on the exon pieces, a thin line across the gaps."""
    ax.plot([segs[0][0], segs[-1][1]], [y, y], color=col, lw=0.8, zorder=2)
    for x0, x1 in segs:
        ax.plot([x0, x1], [y, y], color=col, lw=lw, solid_capstyle="butt", zorder=3)

# (a) isoforms — exons E1..E5, E2 and E4 are optional
a = axes[0]
E = {1: (0.4, 1.3), 2: (1.7, 2.4), 3: (2.8, 3.7), 4: (4.1, 4.8), 5: (5.2, 6.3)}
a.text(0.4, 5.6, "gene", ha="left", va="center", fontsize=9, color=TEXT2)
a.plot([0.4, 6.3], [5.15, 5.15], color=LINE, lw=1)
for k, (x0, x1) in E.items():
    box(a, x0, x1, 5.15, DARK if k in (1, 3, 5) else tint(0.45), label=f"E{k}")
a.text(6.55, 5.15, "E2, E4 optional", ha="left", va="center", fontsize=8.4, color=TEXT2)

a.text(0.4, 4.35, "short reads (150 bp)", ha="left", va="center", fontsize=9, fontweight="semibold", color=TEXT)
pieces = [((1, 2),), ((2, 3),), ((1, 3),), ((3, 4),), ((4, 5),), ((3, 5),)]
for i, ((p, q),) in enumerate(pieces):
    y = 3.9 - (i % 3) * 0.32 - (i // 3) * 0.0
    xo = 0 if i < 3 else 0.0
    read(a, [(E[p][1] - 0.35, E[p][1]), (E[q][0], E[q][0] + 0.35)], y if i < 3 else 3.9 - (i - 3) * 0.32,
         GREYF if True else CRIMSON)
a.text(6.55, 3.6, "each piece covers\n1–2 exons only", ha="left", va="center", fontsize=8.4, color=TEXT2,
       linespacing=1.3)

a.text(0.4, 2.55, "long reads (whole molecule)", ha="left", va="center", fontsize=9, fontweight="semibold",
       color=CRIMSON)
for i, ex in enumerate(((1, 2, 3, 5), (1, 3, 4, 5))):
    read(a, [E[k] for k in ex], 2.05 - i * 0.42, CRIMSON, lw=4)
    a.text(6.55, 2.05 - i * 0.42, "E" + "–E".join(map(str, ex)), ha="left", va="center", fontsize=8.6,
           color=CRIMSON, fontweight="semibold")
a.text(0.4, 0.55, "Short pieces fit E1–E2–E3–E5 + E1–E3–E4–E5 just as well as\nE1–E2–E3–E4–E5 + E1–E3–E5. "
       "A long read shows which exons\nsit on the same molecule.", ha="left", va="center", fontsize=8.6,
       color=TEXT, linespacing=1.35)
a.set_title("(a) Isoforms", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=24)
a.text(0, 1.03, "which exons belong together?", transform=a.transAxes, ha="left", va="bottom", fontsize=9.1,
       color=TEXT2)

# (b) fusion — gene X joined to gene Y
b = axes[1]
X = [(0.6, 1.6), (2.0, 3.0)]; Y = [(3.3, 4.3), (4.7, 5.7)]
b.text(0.6, 5.75, "fusion transcript", ha="left", va="center", fontsize=9, color=TEXT2)
b.text(1.8, 4.75, "gene X", ha="center", va="center", fontsize=8.4, color=DARK, fontweight="semibold")
b.text(4.5, 4.75, "gene Y", ha="center", va="center", fontsize=8.4, color=CRIMSON, fontweight="semibold")
for (x0, x1), lab in zip(X, ("X1", "X2")): box(b, x0, x1, 5.15, DARK, label=lab)
for (x0, x1), lab in zip(Y, ("Y2", "Y3")): box(b, x0, x1, 5.15, tint(0.55), label=lab)
b.plot([3.15, 3.15], [4.85, 5.45], color=CRIMSON, lw=1.2, ls=(0, (2, 2)))
b.text(5.95, 5.15, "dashed = breakpoint", ha="left", va="center", fontsize=8.2, color=CRIMSON)

b.text(0.6, 4.3, "short reads", ha="left", va="center", fontsize=9, fontweight="semibold", color=TEXT)
for i, (x0, x1) in enumerate(((0.7, 1.5), (2.1, 2.9), (3.4, 4.2), (4.8, 5.6))):
    b.plot([x0, x1], [3.9 - (i % 2) * 0.3, 3.9 - (i % 2) * 0.3], color=GREYF, lw=3.2, solid_capstyle="butt")
read(b, [(2.6, 3.0), (3.3, 3.7)], 3.3, CRIMSON)
b.text(3.9, 3.3, "crosses the breakpoint (rare)", ha="left", va="center", fontsize=8.2, color=CRIMSON)
b.plot([2.2, 2.7], [2.92, 2.92], color=CRIMSON, lw=3.2, solid_capstyle="butt")
b.plot([4.9, 5.4], [2.92, 2.92], color=CRIMSON, lw=3.2, solid_capstyle="butt")
b.plot([2.7, 4.9], [2.92, 2.92], color=CRIMSON, lw=0.9, ls=(0, (1.5, 1.5)))
b.text(5.6, 2.92, "paired-end mates in X and Y", ha="left", va="center", fontsize=8.2, color=CRIMSON)
b.text(6.2, 4.05, "most reads look\nlike normal X or Y", ha="left", va="center", fontsize=8.4, color=TEXT2,
       linespacing=1.3)

b.text(0.6, 2.25, "long read", ha="left", va="center", fontsize=9, fontweight="semibold", color=CRIMSON)
read(b, X + Y, 1.8, CRIMSON, lw=4)
b.text(6.0, 1.8, "X–Y in one read", ha="left", va="center", fontsize=8.6, color=CRIMSON, fontweight="semibold")
b.text(0.6, 0.55, "Short reads detect a fusion only from the few reads or pairs that\nstraddle the breakpoint. "
       "A long read shows both partners and\nthe exons kept from each.", ha="left", va="center", fontsize=8.6,
       color=TEXT, linespacing=1.35)
b.set_title("(b) Fusion genes", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=24)
b.text(0, 1.03, "two genes joined in one transcript", transform=b.transAxes, ha="left", va="bottom", fontsize=9.1,
       color=TEXT2)

header(axes[0], "Before you sequence · read length", "Short reads see pieces; a long read sees the whole molecule",
       dy=44, x1=2.08)
fig.savefig(OUT + "long_read.png"); print("ok h28 long read")
