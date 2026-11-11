from gallery import *

# Differential expression vs differential splicing: same total, different isoform mix
fig = plt.figure(figsize=(9.8, 5.6))
gs = fig.add_gridspec(2, 2, width_ratios=[1.6, 1.0], height_ratios=[1, 1], wspace=0.22, hspace=0.75,
                      left=0.01, right=0.98, top=0.8, bottom=0.08)
DARK = "#5E6670"; GREYF = "#B8BEC7"
EX = [(0.6, 1.6), (3.0, 4.0), (5.3, 6.1), (8.1, 9.3)]       # exon 1–4
ISO1, ISO2 = CRIMSON, DARK

def panel(ax, t, sub, pad=22):
    ax.set_title(t, loc="left", fontsize=12.5, fontweight="semibold", color=TEXT, pad=pad)
    ax.text(0, 1.02, sub, transform=ax.transAxes, ha="left", va="bottom", fontsize=9.3, color=TEXT2)

# left: isoforms and sashimi-style coverage ------------------------------------------------------------
L = fig.add_subplot(gs[:, 0]); L.set_xlim(-1.6, 9.6); L.set_ylim(-1.6, 10.9); L.axis("off")
panel(L, "(a) One gene, two isoforms", "exon 3 is included in isoform 1 and skipped in isoform 2")
def iso(y, exons, col, lab):
    L.text(-1.6, y, lab, ha="left", va="center", fontsize=9.3, fontweight="semibold", color=col)
    for (a0, b0), (a1, b1) in zip(exons[:-1], exons[1:]):
        m = (b0 + a1) / 2
        L.plot([b0, m, a1], [y + 0.2, y + 0.55, y + 0.2], color=col, lw=1.1)
    for a0, b0 in exons:
        L.add_patch(mp.Rectangle((a0, y - 0.22), b0 - a0, 0.44, fc=col, ec="none"))
iso(10.1, EX, ISO1, "isoform 1")
iso(8.95, [EX[0], EX[1], EX[3]], ISO2, "isoform 2")
for k, (a0, b0) in enumerate(EX):
    L.text((a0 + b0) / 2, 8.2, f"exon {k + 1}", ha="center", va="center", fontsize=8.5,
           color=CRIMSON if k == 2 else TEXT2, fontweight="semibold" if k == 2 else "normal")

def sashimi(y0, inc, skip, lab):
    L.text(-1.6, y0 + 0.8, lab, ha="left", va="center", fontsize=10, fontweight="semibold", color=TEXT)
    L.text(-1.6, y0 + 0.35, "reads", ha="left", va="center", fontsize=8.5, color=TEXT2)
    L.plot([EX[0][0], EX[3][1]], [y0, y0], color=LINE, lw=0.8)
    h = 1.3
    for k, (a0, b0) in enumerate(EX):
        hh = h * (inc / 100 if k == 2 else 1)
        L.add_patch(mp.Rectangle((a0, y0), b0 - a0, hh, fc=tint(0.5) if k == 2 else GREYF, ec="none"))
    def arc(x0, x1, n, up, col):
        xs = np.linspace(x0, x1, 60); t = (xs - x0) / (x1 - x0)
        bump = (0.55 if up else 0.75) * np.sin(np.pi * t)
        base = y0 + h + 0.12 if up else y0 - 0.1
        L.plot(xs, base + bump if up else base - bump, color=col, lw=0.6 + n / 30, solid_capstyle="round")
        L.text((x0 + x1) / 2, (base + 0.82) if up else (base - 1.02), str(n), ha="center", va="center",
               fontsize=8.8, color=col, fontweight="semibold")
    arc(EX[1][1], EX[2][0], inc, True, ISO1)
    arc(EX[2][1], EX[3][0], inc, True, ISO1)
    arc(EX[1][1], EX[3][0], skip, False, ISO2)
    arc(EX[0][1], EX[1][0], 100, True, GREYF)
sashimi(4.1, 80, 20, "control")
sashimi(0.1, 25, 75, "case")
L.text(-1.6, 7.45, "arcs = reads that jump an intron (number = read count) · bars = reads on each exon",
       ha="left", va="center", fontsize=8.5, color=TEXT2)

# right top: total expression ---------------------------------------------------------------------------
R1 = fig.add_subplot(gs[0, 1])
panel(R1, "(b) Differential expression?", "total reads from the gene")
R1.bar([0, 1], [100, 100], width=0.55, color=[GREYF, tint(0.5)])
for x, v in ((0, 100), (1, 100)):
    R1.text(x, v + 4, f"{v}", ha="center", va="bottom", fontsize=9, color=TEXT)
R1.set_xticks([0, 1]); R1.set_xticklabels(["control", "case"]); R1.set_ylim(0, 135); R1.set_yticks([])
R1.text(1.45, 60, "no change\nlog₂FC = 0", ha="left", va="center", fontsize=10, fontweight="bold", color=DARK,
        linespacing=1.3)

# right bottom: isoform mix ------------------------------------------------------------------------------
R2 = fig.add_subplot(gs[1, 1])
panel(R2, "(c) Differential splicing?", "share of each isoform")
for x, f1 in ((0, 0.8), (1, 0.25)):
    R2.bar(x, f1, width=0.55, color=ISO1); R2.bar(x, 1 - f1, bottom=f1, width=0.55, color=ISO2)
    R2.text(x, f1 / 2, f"{f1:.0%}", ha="center", va="center", fontsize=8.8, color="white", fontweight="semibold")
    R2.text(x, f1 + (1 - f1) / 2, f"{1 - f1:.0%}", ha="center", va="center", fontsize=8.8, color="white",
            fontweight="semibold")
R2.set_xticks([0, 1]); R2.set_xticklabels(["control", "case"]); R2.set_ylim(0, 1); R2.set_yticks([])
R2.text(1.45, 0.62, "yes", ha="left", va="center", fontsize=10.5, fontweight="bold", color=CRIMSON)
R2.text(1.45, 0.32, "exon 3 inclusion\n0.80 → 0.25", ha="left", va="center", fontsize=9, color=CRIMSON,
        linespacing=1.3)
for ax in (R1, R2):
    ax.set_xlim(-0.5, 2.6)
    for sp in ("left",): ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_color(LINE); ax.tick_params(colors=TEXT2, labelsize=9, length=0)

header(L, "Beyond DE · splicing and isoforms",
       "The gene's total can stay flat while its isoform mix flips", dy=40, x1=1.68)
fig.savefig(OUT + "diff_splicing.png"); print("ok h12 splicing")
