from gallery import *

# Gene switches: a TF reads its motif; promoter at the gene, enhancer far away;
# a loop decides which enhancer can reach which promoter.
W, Hh = 11.6, 9.7
fig, ax = plt.subplots(figsize=(11.5, 11.5 * Hh / W + 0.9))
ax.set_xlim(0, W); ax.set_ylim(-1.5, 8.2); ax.set_aspect("equal"); ax.axis("off")
MONO = "DejaVu Sans Mono"; DARK = "#5E6670"; PROT = "#B9BFC8"

def title(x, y, n, t):
    ax.add_patch(mp.Circle((x + 0.17, y), 0.15, fc=CRIMSON, ec="none"))
    ax.text(x + 0.17, y, str(n), ha="center", va="center", fontsize=10, fontweight="bold", color="white")
    ax.text(x + 0.42, y, t, ha="left", va="center", fontsize=12, fontweight="semibold", color=TEXT)

def tf(x, y, w=0.62, h=0.36, fc=CRIMSON, label="TF", fs=9):
    ax.add_patch(mp.FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0,rounding_size=0.15",
                                   fc=fc, ec="none", zorder=5))
    ax.text(x, y, label, ha="center", va="center", fontsize=fs, fontweight="bold", color="white", zorder=6)

# ① TF reads a motif ------------------------------------------------------------------------------
title(0.0, 7.65, 1, "A TF binds only where its motif is")
def strand(y, seq, motif, bound):
    x0, dx = 0.35, 0.29
    if motif:
        a, b = motif
        ax.add_patch(mp.FancyBboxPatch((x0 + a * dx - 0.13, y - 0.17), (b - a) * dx - 0.05, 0.34,
                                       boxstyle="round,pad=0,rounding_size=0.06", fc=tint(0.22), ec="none", zorder=1))
    for i, ch in enumerate(seq):
        inm = motif and motif[0] <= i < motif[1]
        ax.text(x0 + i * dx, y, ch, ha="center", va="center", fontsize=11, family=MONO,
                fontweight="bold" if inm else "normal", color=CRIMSON if inm else TEXT2, zorder=2)
    ax.plot([x0 - 0.15, x0 + (len(seq) - 1) * dx + 0.15], [y - 0.27, y - 0.27], color=LINE, lw=1.6)
    return x0, dx
x0, dx = strand(6.5, "CATGACTCATG", (2, 9), True)
mx = x0 + 5 * dx
tf(mx, 7.03, w=1.6, h=0.4, label="TF  (e.g. AP-1)")
ax.annotate("", (mx, 6.7), (mx, 6.83), arrowprops=dict(arrowstyle="-|>", color=CRIMSON, lw=1.2, mutation_scale=10))
ax.text(x0 + 11 * dx + 0.05, 6.5, "motif TGACTCA\n→ binds", ha="left", va="center", fontsize=9.5,
        color=CRIMSON, fontweight="semibold", linespacing=1.35)
strand(4.95, "CATGCCTGATG", None, False)
tf(mx, 5.72, w=1.6, h=0.4, fc="#C9CDD3", label="TF")
ax.text(mx, 5.3, "✕", ha="center", va="center", fontsize=12, color=DARK)
ax.text(x0 + 11 * dx + 0.05, 4.95, "no motif\n→ no binding", ha="left", va="center", fontsize=9.5,
        color=TEXT2, linespacing=1.35)
ax.text(0.2, 4.3, "TF (transcription factor) = a protein that reads a short DNA word,\n"
        "the motif (~6–12 bp), and switches nearby genes on or off", ha="left", va="center",
        fontsize=9, color=TEXT2, linespacing=1.45)

# ② promoter vs enhancer on the linear DNA ------------------------------------------------------------
X2 = 5.35
ax.plot([X2 - 0.25, X2 - 0.25], [4.0, 7.9], color="#E3E6EB", lw=0.9)
title(X2, 7.65, 2, "Promoter at the gene, enhancer far away")
y = 6.05; ENH = X2 + 0.75; TSS = X2 + 3.35; END = W - 0.15
ax.plot([X2, END], [y, y], color=TEXT2, lw=1.2, zorder=1)
ax.add_patch(mp.Rectangle((ENH - 0.3, y - 0.12), 0.6, 0.24, fc=tint(0.35), ec="none", zorder=2))
for a, b in [(TSS, TSS + 0.32), (TSS + 0.85, TSS + 1.05), (TSS + 1.45, TSS + 1.65), (TSS + 2.05, TSS + 2.25),
             (END - 0.4, END - 0.1)]:
    ax.add_patch(mp.Rectangle((a, y - 0.17), b - a, 0.34, fc=TEXT, ec="none", zorder=2))
ax.plot([TSS, TSS], [y + 0.17, y + 0.5], color=TEXT, lw=1.1)
ax.annotate("", (TSS + 0.5, y + 0.5), (TSS, y + 0.5), arrowprops=dict(arrowstyle="-|>", color=TEXT, lw=1.1,
                                                                    mutation_scale=10))
ax.text(TSS + 0.58, y + 0.5, "TSS", ha="left", va="center", fontsize=8.8, color=TEXT2)
tf(ENH, y + 0.42, w=0.5, h=0.3, fs=8); tf(TSS - 0.38, y + 0.42, w=0.5, h=0.3, fs=8)
ax.add_patch(mp.Rectangle((TSS - 0.75, y - 0.1), 0.75, 0.2, fc=tint(0.6), ec="none", zorder=2))
# labels
ax.text(ENH, y - 0.42, "enhancer", ha="center", va="top", fontsize=10, fontweight="semibold", color=CRIMSON)
ax.text(TSS - 0.38, y - 0.42, "promoter", ha="center", va="top", fontsize=10, fontweight="semibold", color=CRIMSON)
ax.text((TSS + END) / 2 + 0.2, y - 0.42, "gene body", ha="center", va="top", fontsize=10, fontweight="semibold",
        color=TEXT)
# distance bracket
yb = y + 0.98
ax.annotate("", (ENH, yb), (TSS - 0.38, yb), arrowprops=dict(arrowstyle="<->", color=TEXT2, lw=1))
ax.text((ENH + TSS - 0.38) / 2, yb + 0.14, "kb – Mb apart", ha="center", va="bottom", fontsize=9, color=TEXT2)
ax.text(X2, 4.55, "promoter = the switch right at the start of the gene (around the TSS)\n"
        "enhancer = a remote booster: can sit far away, upstream or downstream,\n"
        "even inside another gene · TFs bind both", ha="left", va="center", fontsize=9, color=TEXT2,
        linespacing=1.45)

ax.plot([0, W], [3.75, 3.75], color="#E3E6EB", lw=0.9)

# ③ the loop decides who can reach whom --------------------------------------------------------------
title(0.0, 3.35, 3, "In 3D, a loop decides which enhancer can reach which gene")
cx, cy, R, y0 = 3.55, 1.05, 1.3, -0.85        # loop circle and baseline
A0, SPAN = 250, 320                              # loop starts at 250°, runs clockwise 320°
LOOP_KB = 125
def ang(kb): return np.deg2rad(A0 - kb / LOOP_KB * SPAN)
def pt(kb): return cx + R * np.cos(ang(kb)), cy + R * np.sin(ang(kb))
th = np.linspace(ang(0), ang(LOOP_KB), 400)
ax.plot(cx + R * np.cos(th), cy + R * np.sin(th), color=LINE, lw=3, zorder=1, solid_capstyle="round")
pL, pR = pt(0), pt(LOOP_KB)
nL, nR = (cx - 0.13, y0), (cx + 0.13, y0)
ax.plot([pL[0], nL[0]], [pL[1], nL[1]], color=LINE, lw=3, zorder=1)
ax.plot([pR[0], nR[0]], [pR[1], nR[1]], color=LINE, lw=3, zorder=1)
ax.plot([0.25, nL[0]], [y0, y0], color=LINE, lw=3, zorder=1, solid_capstyle="round")
ax.plot([nR[0], 6.9], [y0, y0], color=LINE, lw=3, zorder=1, solid_capstyle="round")
# cohesin ring + CTCF at the neck
ax.add_patch(mp.Ellipse((cx, y0 + 0.42), 0.75, 0.32, fc="none", ec=DARK, lw=2.2, zorder=4))
for xx in (cx - 0.3, cx + 0.3):
    ax.add_patch(mp.Circle((xx, y0 + 0.12), 0.1, fc=PROT, ec="none", zorder=3))
ax.text(cx + 0.5, y0 + 0.62, "cohesin", ha="left", va="center", fontsize=8.5, color=DARK)
ax.text(cx + 0.5, y0 + 0.3, "CTCF (loop anchors)", ha="left", va="center", fontsize=8.5, color=DARK)
# gene inside the loop: promoter at 15 kb, body 15–45 kb
gb = np.linspace(ang(15), ang(45), 60)
ax.plot(cx + R * np.cos(gb), cy + R * np.sin(gb), color=TEXT, lw=5, zorder=2, solid_capstyle="butt")
P = pt(15); E1 = pt(110)
ax.add_patch(mp.Circle(P, 0.13, fc=CRIMSON, ec="white", lw=1.2, zorder=4))
ax.text(P[0] - 0.22, P[1] - 0.05, "promoter", ha="right", va="center", fontsize=9.5, color=CRIMSON,
        fontweight="semibold")
G = pt(36)
ax.text(G[0] - 0.22, G[1], "gene", ha="right", va="center", fontsize=9.5, color=TEXT, fontweight="semibold")
ax.add_patch(mp.Circle(E1, 0.16, fc=tint(0.45), ec="white", lw=1.2, zorder=4))
ax.text(E1[0] + 0.25, E1[1] + 0.05, "enhancer 1", ha="left", va="center", fontsize=9.5, color=CRIMSON,
        fontweight="semibold")
ax.annotate("", (P[0] + 0.15, P[1] + 0.05), (E1[0] - 0.18, E1[1] + 0.02),
            arrowprops=dict(arrowstyle="-|>", color=CRIMSON, lw=1.6, mutation_scale=13,
                            connectionstyle="arc3,rad=0.25"))
ax.text(cx + 0.1, cy - 0.1, "✓ same loop\n→ they meet", ha="center", va="center", fontsize=9.5,
        color=CRIMSON, fontweight="semibold", linespacing=1.3)
# enhancer 2: outside the loop, close on the line
E2 = (1.25, y0)
ax.add_patch(mp.Circle(E2, 0.16, fc="#C9CDD3", ec="white", lw=1.2, zorder=4))
ax.text(E2[0], y0 - 0.32, "enhancer 2", ha="center", va="top", fontsize=9.5, color=DARK, fontweight="semibold")
tip = (P[0] - 0.55, P[1] - 0.45)
ax.annotate("", tip, (E2[0] + 0.1, E2[1] + 0.16),
            arrowprops=dict(arrowstyle="-", color=DARK, lw=1.3, ls=(0, (3, 2)),
                            connectionstyle="arc3,rad=-0.2"))
ax.text(tip[0] + 0.04, tip[1] + 0.06, "✕", ha="center", va="center", fontsize=16, color=DARK, fontweight="bold")
ax.text(0.15, 0.75, "outside the loop\n→ cannot reach", ha="left", va="center", fontsize=9.5, color=DARK,
        fontweight="semibold", linespacing=1.3)

# the same locus on the linear genome
LX0, LX1 = 7.35, W - 0.2
KB0, KB1 = -35, 135
def lx(kb): return LX0 + (kb - KB0) / (KB1 - KB0) * (LX1 - LX0)
ax.text(LX0, 2.6, "The same stretch laid out straight", ha="left", va="center", fontsize=10.5,
        fontweight="semibold", color=TEXT)
yl = 1.0
ax.plot([lx(KB0), lx(KB1)], [yl, yl], color=TEXT2, lw=1.2, zorder=1)
ax.add_patch(mp.Rectangle((lx(15), yl - 0.13), lx(45) - lx(15), 0.26, fc=TEXT, ec="none", zorder=2))
ax.add_patch(mp.Circle((lx(15), yl), 0.1, fc=CRIMSON, ec="white", lw=1, zorder=3))
ax.add_patch(mp.Circle((lx(110), yl), 0.12, fc=tint(0.45), ec="white", lw=1, zorder=3))
ax.add_patch(mp.Circle((lx(-15), yl), 0.12, fc="#C9CDD3", ec="white", lw=1, zorder=3))
for kb, d in ((0, 1), (LOOP_KB, -1)):            # convergent CTCF sites
    x = lx(kb)
    ax.add_patch(mp.Polygon([(x - 0.09 * d, yl - 0.13), (x - 0.09 * d, yl + 0.13), (x + 0.11 * d, yl)],
                            closed=True, fc=DARK, ec="none", zorder=3))
# loop arc over the anchors
xs = np.linspace(lx(0), lx(LOOP_KB), 100)
ax.plot(xs, yl + 0.25 + 0.6 * np.sin(np.pi * (xs - xs[0]) / (xs[-1] - xs[0])), color=LINE, lw=1.6)
ax.text((lx(0) + lx(LOOP_KB)) / 2, yl + 0.92, "one loop", ha="center", va="bottom", fontsize=9, color=TEXT2)
ax.text(lx(-15), yl - 0.25, "E2", ha="center", va="top", fontsize=9, color=DARK, fontweight="semibold")
ax.text(lx(15), yl - 0.25, "P", ha="center", va="top", fontsize=9, color=CRIMSON, fontweight="semibold")
ax.text(lx(30), yl + 0.2, "gene", ha="center", va="bottom", fontsize=8.5, color=TEXT)
ax.text(lx(110), yl - 0.25, "E1", ha="center", va="top", fontsize=9, color=CRIMSON, fontweight="semibold")
# distances
def dist(a, b, yy, lab, c, w="normal"):
    ax.annotate("", (lx(a), yy), (lx(b), yy), arrowprops=dict(arrowstyle="<->", color=c, lw=1))
    ax.text((lx(a) + lx(b)) / 2, yy - 0.12, lab, ha="center", va="top", fontsize=9, color=c, fontweight=w)
dist(-15, 15, yl - 0.75, "30 kb\nnear, but ✕", DARK, "semibold")
dist(15, 110, yl - 0.75, "95 kb\nfar, but ✓", CRIMSON, "semibold")
ax.text(LX0, -1.0, "Distance on the line does not decide it.\nBeing inside the same loop does.",
        ha="left", va="center", fontsize=9.5, color=TEXT, linespacing=1.45)

header(ax, "Gene switches · promoter, enhancer, transcription factor",
       "TFs flip the switches; loops decide which enhancer reaches which gene", dy=0)
fig.savefig(OUT + "gene_switch.png"); print("ok g23 gene switch")
