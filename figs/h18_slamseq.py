from gallery import *

# SLAM-seq: 4sU marks new RNA; chemistry turns it into a T>C mismatch in the read; new/old over time → half-life
W, H = 9.8, 4.3
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, (H - 1.0) / H]); ax.set_xlim(0, W); ax.set_ylim(0, H - 1.0); ax.axis("off")
BG = "#F8F9FA"; ARROW = "#9AA1AB"; DARK = "#5E6670"; GREYF = "#C3C8D0"
steps = [("4sU pulse", "new RNA takes up 4sU in place of some U"),
         ("4sU reads as C", "one library, no pull-down"),
         ("New vs old, over time", "→ synthesis and half-life per gene")]
CW, GAP, TOP, CH = 3.08, 0.2, H - 1.12, 2.45
X0 = [0.08 + k * (CW + GAP) for k in range(3)]
for k, (name, sub) in enumerate(steps):
    x = X0[k]
    ax.add_patch(mp.FancyBboxPatch((x, TOP - CH), CW, CH, boxstyle="round,pad=0,rounding_size=0.12",
                                   fc=BG, ec=LINE, lw=0.9))
    ax.scatter([x + 0.27], [TOP - 0.28], s=300, color=tint(0.16), lw=0, zorder=3)
    ax.text(x + 0.27, TOP - 0.28, str(k + 1), ha="center", va="center", fontsize=9.5, fontweight="bold",
            color=CRIMSON, zorder=4)
    ax.text(x + 0.5, TOP - 0.28, name, ha="left", va="center", fontsize=11, fontweight="bold", color=TEXT)
    ax.text(x + 0.16, TOP - 0.6, sub, ha="left", va="center", fontsize=8.6, color=TEXT2)
    if k < 2:
        ax.annotate("", (X0[k + 1] - 0.02, TOP - CH / 2), (x + CW + 0.02, TOP - CH / 2),
                    arrowprops=dict(arrowstyle="-|>", color=ARROW, lw=1.1, mutation_scale=9))

# 1 — new RNA carries 4sU marks; old RNA does not
x = X0[0]
for i, (label, marks) in enumerate((("new RNA", (0.55, 1.15, 1.75)), ("new RNA", (0.85, 1.55)),
                                    ("old RNA", ()), ("old RNA", ()))):
    yy = TOP - 1.0 - i * 0.33
    ax.plot([x + 0.25, x + 2.15], [yy, yy], color=DARK, lw=2.6, solid_capstyle="butt")
    for m in marks:
        ax.plot([x + 0.25 + m, x + 0.25 + m], [yy - 0.1, yy + 0.1], color=CRIMSON, lw=2.4)
    ax.text(x + 2.3, yy, label, ha="left", va="center", fontsize=8.2,
            color=CRIMSON if marks else TEXT2, fontweight="semibold" if marks else "normal")
ax.text(x + 0.2, TOP - CH + 0.2, "crimson tick = 4sU", ha="left", va="center", fontsize=8.3, color=TEXT2)

# 2 — reference vs reads: the new read shows a T>C mismatch
x = X0[1]
ref = "ACTGATCAGT"; dx = 0.2; xs0 = x + 0.45
rows = (("genome", ref, None), ("old read", ref, None), ("new read", ref[:5] + "C" + ref[6:], 5))
for r, (lab, seq, hit) in enumerate(rows):
    yy = TOP - 1.05 - r * 0.42
    ax.text(x + 0.2, yy, lab, ha="left", va="center", fontsize=8.2, color=TEXT2 if r != 2 else CRIMSON,
            fontweight="semibold" if r == 2 else "normal")
    for j, b in enumerate(seq):
        bx = xs0 + 0.55 + j * dx
        if hit is not None and j == hit:
            ax.add_patch(mp.Rectangle((bx - 0.1, yy - 0.14), 0.2, 0.28, fc=CRIMSON, ec="none"))
        ax.text(bx, yy, b, ha="center", va="center", fontsize=9.5, family="DejaVu Sans Mono",
                color="white" if (hit is not None and j == hit) else (TEXT if r else TEXT2),
                fontweight="bold" if (hit is not None and j == hit) else "normal")
ax.text(x + 0.2, TOP - CH + 0.2, "T in the genome, C in the read → new RNA", ha="left", va="center",
        fontsize=8.3, color=TEXT)

# 3 — fraction of reads with T>C rises fast for an unstable RNA, slowly for a stable one
x = X0[2]; px0, px1, py0, py1 = x + 0.45, x + 2.85, TOP - 2.1, TOP - 0.95
ax.plot([px0, px1], [py0, py0], color=LINE, lw=1); ax.plot([px0, px0], [py0, py1], color=LINE, lw=1)
t = np.linspace(0, 1, 80)
for k, col, lab in ((7.0, CRIMSON, "short half-life"), (1.2, GREYF, "long half-life")):
    f = 1 - np.exp(-k * t)
    ax.plot(px0 + t * (px1 - px0), py0 + f * (py1 - py0) * 0.92, color=col, lw=2.2)
    ax.text(px0 + 0.62 * (px1 - px0) + 0.05, py0 + (1 - np.exp(-k * 0.62)) * (py1 - py0) * 0.92 + (0.1 if k > 2 else -0.14),
            lab, ha="left", va="center", fontsize=8.2, color=CRIMSON if k > 2 else TEXT2,
            fontweight="semibold" if k > 2 else "normal")
ax.text(px0 - 0.08, (py0 + py1) / 2, "% new", ha="right", va="center", fontsize=8, color=TEXT2, rotation=90)
ax.text(px1, py0 - 0.14, "hours of labelling", ha="right", va="center", fontsize=8, color=TEXT2)

bx0, bx1, by = 0.08, W - 0.08, 0.12
ax.add_patch(mp.FancyBboxPatch((bx0, by), bx1 - bx0, 0.46, boxstyle="round,pad=0,rounding_size=0.08",
                               fc=CRIMSON, ec="none"))
ax.text((bx0 + bx1) / 2, by + 0.23, "T>C marks split each gene's reads into new and old RNA → synthesis and decay, separately",
        ha="center", va="center", fontsize=10.3, fontweight="bold", color="white")
header(ax, "Beyond RNA-seq · metabolic labelling (SLAM-seq)", "Tell new RNA from old RNA inside one library")
fig.savefig(OUT + "slam_seq.png"); print("ok h18 slam-seq")
