from gallery import *

# How bulk RNA-seq works: remove rRNA (poly(A) capture), make cDNA, sequence, map and count
W, H = 9.8, 4.5
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, (H - 1.0) / H]); ax.set_xlim(0, W); ax.set_ylim(0, H - 1.0); ax.axis("off")
BG = "#F8F9FA"; ARROW = "#9AA1AB"; DARK = "#5E6670"; GREYF = "#C3C8D0"
rng = np.random.default_rng(3)

steps = [("RNA in the tube", "rRNA ≈ 80–90% of it"),
         ("Keep poly(A) RNA", "oligo-dT beads"),
         ("cDNA library", "fragment, copy, add adapters"),
         ("Sequence", "millions of short reads"),
         ("Map and count", "reads per gene")]
CW, GAP, TOP, CH = 1.8, 0.17, H - 1.12, 2.62
X0 = [0.08 + k * (CW + GAP) for k in range(5)]
for k, (name, sub) in enumerate(steps):
    x = X0[k]
    ax.add_patch(mp.FancyBboxPatch((x, TOP - CH), CW, CH, boxstyle="round,pad=0,rounding_size=0.12",
                                   fc=BG, ec=LINE, lw=0.9))
    ax.scatter([x + 0.25], [TOP - 0.27], s=280, color=tint(0.16), lw=0, zorder=3)
    ax.text(x + 0.25, TOP - 0.27, str(k + 1), ha="center", va="center", fontsize=9.5, fontweight="bold",
            color=CRIMSON, zorder=4)
    ax.text(x + 0.47, TOP - 0.27, name, ha="left", va="center", fontsize=10.3, fontweight="bold", color=TEXT)
    ax.text(x + 0.14, TOP - 0.58, sub, ha="left", va="center", fontsize=8.3, color=TEXT2)
    if k < 4:
        ax.annotate("", (X0[k + 1] - 0.02, TOP - CH / 2), (x + CW + 0.02, TOP - CH / 2),
                    arrowprops=dict(arrowstyle="-|>", color=ARROW, lw=1.1, mutation_scale=9))

def mrna(x, y, L=0.95, lw=2.6, tail=True):
    ax.plot([x, x + L], [y, y], color=DARK, lw=lw, solid_capstyle="butt", zorder=4)
    if tail:
        ax.text(x + L + 0.03, y, "AAAA", ha="left", va="center", fontsize=6.6, color=CRIMSON,
                fontweight="bold", zorder=4)

# 1 — the tube: lots of rRNA, a few mRNAs with tails
x = X0[0]
for i in range(7):
    xx = x + 0.15 + rng.uniform(0, 0.55); yy = TOP - 0.88 - i * 0.13
    ax.plot([xx, xx + rng.uniform(0.7, 0.95)], [yy, yy], color=GREYF, lw=2.4, solid_capstyle="butt")
mrna(x + 0.15, TOP - 1.8, L=0.85); mrna(x + 0.3, TOP - 1.96, L=0.7)
ax.text(x + 0.14, TOP - CH + 0.42, "grey = rRNA", ha="left", va="center", fontsize=8, color=TEXT2)
ax.text(x + 0.14, TOP - CH + 0.2, "dark = mRNA + poly(A) tail", ha="left", va="center", fontsize=8, color=TEXT)

# 2 — oligo-dT bead captures tails; rRNA washed away
x = X0[1]; cy = TOP - 1.45
ax.add_patch(mp.Circle((x + 1.42, cy), 0.25, fc=tint(0.2), ec=CRIMSON, lw=1.1, zorder=3))
ax.text(x + 1.42, cy, "bead", ha="center", va="center", fontsize=7.4, color=CRIMSON, zorder=4)
for dy in (0.18, -0.18):
    ax.plot([x + 0.12, x + 0.82], [cy + dy, cy + dy], color=DARK, lw=2.6, solid_capstyle="butt", zorder=4)
    ax.text(x + 0.84, cy + dy, "AAAA", ha="left", va="center", fontsize=6.6, color=CRIMSON, fontweight="bold")
    ax.plot([x + 1.08, x + 1.18], [cy + dy, cy + dy], color=CRIMSON, lw=0.8)
ax.text(x + 0.14, TOP - CH + 0.62, "TTTT on the bead pairs", ha="left", va="center", fontsize=8, color=TEXT)
ax.text(x + 0.14, TOP - CH + 0.41, "with AAAA on the mRNA", ha="left", va="center", fontsize=8, color=TEXT)
ax.text(x + 0.14, TOP - CH + 0.2, "no tail (rRNA) → washed away", ha="left", va="center", fontsize=8,
        color=TEXT2)

# 3 — fragments → cDNA with adapters at both ends
x = X0[2]
for i in range(5):
    yy = TOP - 1.0 - i * 0.24; xx = x + 0.28 + (i % 2) * 0.22
    ax.plot([xx, xx + 0.12], [yy, yy], color=CRIMSON, lw=3, solid_capstyle="butt")
    ax.plot([xx + 0.12, xx + 0.95], [yy, yy], color=DARK, lw=3, solid_capstyle="butt")
    ax.plot([xx + 0.95, xx + 1.07], [yy, yy], color=CRIMSON, lw=3, solid_capstyle="butt")
ax.text(x + 0.14, TOP - CH + 0.41, "RNA copied into cDNA", ha="left", va="center", fontsize=8,
        color=TEXT)
ax.text(x + 0.14, TOP - CH + 0.2, "crimson ends = adapters", ha="left", va="center", fontsize=8, color=TEXT2)

# 4 — reads
x = X0[3]
for i in range(9):
    yy = TOP - 0.92 - i * 0.135; xx = x + 0.18 + rng.uniform(0, 0.75)
    ax.plot([xx, xx + 0.5], [yy, yy], color=tint(0.6), lw=2.2, solid_capstyle="butt")
ax.text(x + 0.14, TOP - CH + 0.41, "each read = the sequence of", ha="left", va="center", fontsize=8, color=TEXT)
ax.text(x + 0.14, TOP - CH + 0.2, "one end of one fragment", ha="left", va="center", fontsize=8, color=TEXT)

# 5 — map to genes, count
x = X0[4]
for g, (yy, n, cnt) in enumerate(((TOP - 1.5, 7, "120"), (TOP - 2.0, 2, "15"))):
    ax.plot([x + 0.12, x + 1.25], [yy, yy], color=LINE, lw=1)
    ax.add_patch(mp.Rectangle((x + 0.25, yy - 0.07), 0.9, 0.14, fc=DARK, ec="none", zorder=3))
    for i in range(n):
        xx = x + 0.25 + (i * 0.37) % 0.6; ry = yy + 0.13 + (i // 2) * 0.075
        ax.plot([xx, xx + 0.28], [ry, ry], color=tint(0.6), lw=1.8, solid_capstyle="butt")
    ax.text(x + 0.12, yy - 0.2, "gene A" if g == 0 else "gene B", ha="left", va="center", fontsize=8,
            color=TEXT2)
    ax.text(x + 1.68, yy, cnt, ha="right", va="center", fontsize=10.5, fontweight="bold",
            color=CRIMSON if g == 0 else TEXT)
ax.text(x + 1.68, TOP - 1.08, "count", ha="right", va="center", fontsize=8, color=TEXT2)
ax.text(x + 0.14, TOP - CH + 0.2, "→ genes × samples table", ha="left", va="center", fontsize=8, color=TEXT)

bx0, bx1, by = 0.08, W - 0.08, 0.12
ax.add_patch(mp.FancyBboxPatch((bx0, by), bx1 - bx0, 0.46, boxstyle="round,pad=0,rounding_size=0.08",
                               fc=CRIMSON, ec="none"))
ax.text((bx0 + bx1) / 2, by + 0.23, "More copies of an RNA in the tube → more fragments → more reads → a higher count",
        ha="center", va="center", fontsize=10.3, fontweight="bold", color="white")
header(ax, "Bulk RNA-seq · how it works", "RNA-seq turns RNA into reads, then counts reads per gene")
fig.savefig(OUT + "rnaseq_principle.png"); print("ok h16 principle")
