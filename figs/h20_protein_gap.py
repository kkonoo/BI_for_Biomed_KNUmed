from gallery import *

# Three steps between mRNA and working protein that RNA-seq cannot see
W, H = 9.8, 4.5
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, (H - 1.0) / H]); ax.set_xlim(0, W); ax.set_ylim(0, H - 1.0); ax.axis("off")
BG = "#F8F9FA"; DARK = "#5E6670"; GREYF = "#C3C8D0"
cards = [("Translation efficiency", "proteins made per mRNA"),
         ("Protein half-life", "minutes for some, weeks for others"),
         ("Post-translational control", "modification, cleavage, localisation")]
CW, GAP, TOP, CH = 3.08, 0.2, H - 1.12, 2.62
X0 = [0.08 + k * (CW + GAP) for k in range(3)]
for k, (name, sub) in enumerate(cards):
    x = X0[k]
    ax.add_patch(mp.FancyBboxPatch((x, TOP - CH), CW, CH, boxstyle="round,pad=0,rounding_size=0.12",
                                   fc=BG, ec=LINE, lw=0.9))
    ax.text(x + 0.16, TOP - 0.28, name, ha="left", va="center", fontsize=11, fontweight="bold", color=TEXT)
    ax.text(x + 0.16, TOP - 0.58, sub, ha="left", va="center", fontsize=8.6, color=TEXT2)

def ribosome(cx, y):
    ax.add_patch(mp.Ellipse((cx, y + 0.08), 0.24, 0.15, fc=DARK, ec="white", lw=0.7, zorder=4))
    ax.add_patch(mp.Ellipse((cx, y + 0.21), 0.17, 0.11, fc=DARK, ec="white", lw=0.7, zorder=4))

def proteins(x, y, n, col=CRIMSON, per_row=6):
    for i in range(n):
        ax.add_patch(mp.Circle((x + (i % per_row) * 0.2, y - (i // per_row) * 0.2), 0.075, fc=col, ec="none",
                               zorder=4))

# 1 — same mRNA, different ribosome load
x = X0[0]
for r, (lab, nrib, nprot) in enumerate((("mRNA P", 4, 9), ("mRNA Q", 1, 2))):
    yy = TOP - 1.15 - r * 0.85
    ax.plot([x + 0.2, x + 1.45], [yy, yy], color=DARK, lw=2.2, solid_capstyle="butt")
    for i in range(nrib): ribosome(x + 0.4 + i * 0.28, yy)
    ax.text(x + 0.2, yy - 0.18, lab, ha="left", va="center", fontsize=8, color=TEXT2)
    ax.annotate("", (x + 1.78, yy + 0.05), (x + 1.52, yy + 0.05),
                arrowprops=dict(arrowstyle="-|>", color="#9AA1AB", lw=1, mutation_scale=8))
    proteins(x + 1.95, yy + 0.15, nprot, col=CRIMSON if r == 0 else tint(0.45), per_row=5)
ax.text(x + 0.16, TOP - CH + 0.22, "same mRNA level, more ribosomes → more protein", ha="left", va="center",
        fontsize=8.2, color=TEXT)

# 2 — same synthesis, different decay
x = X0[1]
for r, (lab, nprot, gone) in enumerate((("stable", 10, 0), ("unstable", 2, 4))):
    yy = TOP - 1.05 - r * 0.85
    ax.text(x + 0.2, yy - 0.05, lab, ha="left", va="center", fontsize=8.6, fontweight="semibold",
            color=TEXT)
    proteins(x + 1.05, yy + 0.05, nprot, col=CRIMSON if r == 0 else tint(0.45), per_row=5)
    for i in range(gone):
        gx = x + 1.6 + i * 0.28
        ax.plot([gx - 0.05, gx + 0.05], [yy - 0.04, yy + 0.06], color=GREYF, lw=1.8)
        ax.plot([gx + 0.02, gx + 0.1], [yy - 0.08, yy - 0.02], color=GREYF, lw=1.8)
    if gone:
        ax.text(x + 1.55, yy - 0.3, "degraded", ha="left", va="center", fontsize=8, color=TEXT2)
ax.text(x + 0.16, TOP - CH + 0.22, "same production, faster decay → less protein", ha="left", va="center",
        fontsize=8.2, color=TEXT)

# 3 — same amount of protein, different activity
x = X0[2]
for r, (lab, nmod) in enumerate((("before", 0), ("after signal", 3))):
    yy = TOP - 1.05 - r * 0.62
    ax.text(x + 0.2, yy, lab, ha="left", va="center", fontsize=8.6, fontweight="semibold", color=TEXT)
    for i in range(4):
        cx = x + 1.45 + i * 0.36
        on = i < nmod
        ax.add_patch(mp.Circle((cx, yy), 0.11, fc=CRIMSON if on else GREYF, ec="none", zorder=4))
        if on:
            ax.add_patch(mp.Circle((cx + 0.1, yy + 0.11), 0.065, fc="white", ec=CRIMSON, lw=0.9, zorder=5))
            ax.text(cx + 0.1, yy + 0.11, "P", ha="center", va="center", fontsize=5.6, color=CRIMSON,
                    fontweight="bold", zorder=6)
ax.text(x + 0.2, TOP - 2.08, "P = phosphorylation (switched on)", ha="left", va="center", fontsize=8, color=TEXT2)
ax.text(x + 0.16, TOP - CH + 0.22, "same amount, different activity", ha="left", va="center", fontsize=8.2,
        color=TEXT)

bx0, bx1, by = 0.08, W - 0.08, 0.12
ax.add_patch(mp.FancyBboxPatch((bx0, by), bx1 - bx0, 0.46, boxstyle="round,pad=0,rounding_size=0.08",
                               fc=CRIMSON, ec="none"))
ax.text((bx0 + bx1) / 2, by + 0.23, "All three happen after the mRNA → RNA-seq cannot see them; measure protein directly (W6)",
        ha="center", va="center", fontsize=10.3, fontweight="bold", color="white")
header(ax, "Bulk RNA-seq · what it cannot see", "Between mRNA and working protein, three more steps")
fig.savefig(OUT + "mrna_to_protein.png"); print("ok h20 protein gap")
