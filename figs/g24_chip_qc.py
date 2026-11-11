from gallery import *
rng = np.random.default_rng(24)

# ChIP-seq QC on one page: FRiP, strand cross-correlation, library complexity, ENCODE blacklist
fig = plt.figure(figsize=(9.8, 9.0))
gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.12, left=0.01, right=0.99, top=0.86, bottom=0.03)
DARK = "#5E6670"; GREYF = "#B8BEC7"

def panel(ax, t, sub):
    ax.set_title(t, loc="left", fontsize=12.5, fontweight="semibold", color=TEXT, pad=22)
    ax.text(0, 1.02, sub, transform=ax.transAxes, ha="left", va="bottom", fontsize=9.3, color=TEXT2)

# (a) FRiP -------------------------------------------------------------------------------------------
a = fig.add_subplot(gs[0, 0]); a.set_xlim(0, 10); a.set_ylim(-0.2, 4.2); a.axis("off")
panel(a, "(a) FRiP", "reads in peaks (shaded) / all reads")
peaks = [(1.3, 2.1), (4.6, 5.3), (7.6, 8.5)]
def reads(y, n_in, n_out, lab, frip, strong):
    for p0, p1 in peaks:
        a.add_patch(mp.Rectangle((p0, y - 0.05), p1 - p0, 1.0, fc=FAINT, ec="none", zorder=0))
    xs_in = np.concatenate([rng.uniform(p0 + 0.05, p1 - 0.05, n_in // 3) for p0, p1 in peaks])
    xs_out = rng.uniform(0.1, 9.9, n_out)
    xs_out = xs_out[~np.any([(xs_out > p0) & (xs_out < p1) for p0, p1 in peaks], axis=0)]
    for xx in xs_out:
        h = rng.uniform(0.05, 0.62); a.plot([xx, xx], [y + h, y + h + 0.16], color=GREYF, lw=1.3)
    for xx in xs_in:
        h = rng.uniform(0.05, 0.75); a.plot([xx, xx], [y + h, y + h + 0.16], color=CRIMSON, lw=1.3)
    a.plot([0, 10], [y - 0.05, y - 0.05], color=LINE, lw=0.9)
    a.text(0, y + 1.12, lab, ha="left", va="bottom", fontsize=9.8, fontweight="semibold",
           color=TEXT if strong else TEXT2)
    a.text(10, y + 1.12, f"FRiP ≈ {frip}", ha="right", va="bottom", fontsize=10, fontweight="semibold",
           color=CRIMSON if strong else TEXT2)
reads(2.55, 60, 45, "strong IP: reads gather in peaks", "55%", True)
reads(0.35, 6, 90, "weak IP: reads scattered like background", "6%", False)
a.text(5, -0.35, "low FRiP → weak antibody or failed IP", ha="center", va="top",
       fontsize=9.2, color=TEXT2)

# (b) strand cross-correlation ------------------------------------------------------------------------
b = fig.add_subplot(gs[0, 1]); b.axis("off")
panel(b, "(b) Strand cross-correlation", "slide + strand reads onto − strand reads")
c1 = b.inset_axes([0, 0.58, 1, 0.42]); c1.set_xlim(0, 10); c1.set_ylim(-0.3, 2.0); c1.axis("off")
c1.plot([0.3, 9.7], [0.2, 0.2], color=LINE, lw=1)
c1.add_patch(mp.FancyBboxPatch((4.2, 0.06), 1.6, 0.3, boxstyle="round,pad=0,rounding_size=0.12",
                               fc=CRIMSON, ec="none", zorder=3))
c1.text(5, 0.21, "bound site", ha="center", va="center", fontsize=8, color="white", fontweight="bold", zorder=4)
for k in range(7):
    xp = rng.normal(3.4, 0.45); yy = 0.6 + k * 0.17
    c1.annotate("", (xp + 0.55, yy), (xp, yy), arrowprops=dict(arrowstyle="-|>", color=CRIMSON, lw=1.1,
                                                              mutation_scale=7))
    xm = rng.normal(6.6, 0.45)
    c1.annotate("", (xm - 0.55, yy), (xm, yy), arrowprops=dict(arrowstyle="-|>", color=DARK, lw=1.1,
                                                              mutation_scale=7))
c1.text(3.4, 1.95, "+ strand reads", ha="center", va="top", fontsize=8.8, color=CRIMSON)
c1.text(6.6, 1.95, "− strand reads", ha="center", va="top", fontsize=8.8, color=DARK)
c1.annotate("", (6.6, -0.12), (3.4, -0.12), arrowprops=dict(arrowstyle="<->", color=TEXT2, lw=0.9))
c1.text(5, -0.2, "≈ fragment length", ha="center", va="top", fontsize=8.8, color=TEXT2)
c2 = b.inset_axes([0.08, 0.0, 0.9, 0.44])
s = np.linspace(0, 500, 500)
good = 0.30 + 0.20 * np.exp(-0.5 * ((s - 200) / 45) ** 2) + 0.07 * np.exp(-0.5 * ((s - 40) / 6) ** 2) \
       - 0.00012 * s
bad = 0.30 + 0.025 * np.exp(-0.5 * ((s - 200) / 50) ** 2) + 0.09 * np.exp(-0.5 * ((s - 40) / 6) ** 2) \
      - 0.00008 * s
c2.plot(s, good, color=CRIMSON, lw=1.8); c2.plot(s, bad, color=TEXT2, lw=1.5)
c2.text(200, good.max() + 0.012, "good: peak at fragment length", ha="center", va="bottom", fontsize=8.8,
        color=CRIMSON, fontweight="semibold")
c2.text(50, 0.43, "phantom peak\n(read length)", ha="left", va="center", fontsize=8.5, color=TEXT2,
        linespacing=1.2)
c2.text(330, bad[330] + 0.012, "failed: phantom\npeak only", ha="left", va="bottom", fontsize=8.8,
        color=TEXT2)
c2.set_xlim(0, 500); c2.set_ylim(0.22, 0.56); c2.set_yticks([])
c2.set_xlabel("strand shift (bp)", color=TEXT2, fontsize=9)
c2.set_ylabel("correlation", color=TEXT2, fontsize=9)
for sp in ("left", "bottom"): c2.spines[sp].set_color(LINE)
c2.tick_params(colors=TEXT2, labelsize=8.5)
c2.text(500, 0.555, "NSC, RSC =\nhow much the\nfragment peak\nstands out", ha="right", va="top", fontsize=8.8,
        color=TEXT, linespacing=1.25)

# (c) library complexity -----------------------------------------------------------------------------
c = fig.add_subplot(gs[1, 0]); c.set_xlim(0, 10); c.set_ylim(-0.6, 4.2); c.axis("off")
panel(c, "(c) Library complexity", "how many different fragments vs PCR copies of the same few")
def library(y, groups, lab, txt, good):
    c.plot([0, 10], [y, y], color=LINE, lw=0.9)
    for x0, n in groups:
        for k in range(n):
            col = CRIMSON if k == 0 else (tint(0.3) if good else GREYF)
            c.add_patch(mp.Rectangle((x0, y + 0.08 + k * 0.15), 0.62, 0.1, fc=col, ec="none"))
    c.text(0, y + 1.3, lab, ha="left", va="bottom", fontsize=9.8, fontweight="semibold",
           color=TEXT if good else TEXT2)
    c.text(10, y + 1.3, txt, ha="right", va="bottom", fontsize=10, fontweight="semibold",
           color=CRIMSON if good else TEXT2)
library(2.35, [(x, 1) for x in np.linspace(0.2, 9.1, 12)], "high: 12 reads = 12 different fragments",
        "NRF = 12/12", True)
library(0.0, [(1.3, 7), (4.6, 3), (7.8, 2)], "low: 12 reads = 3 fragments, copied by PCR",
        "NRF = 3/12", False)
c.add_patch(mp.Rectangle((0.2, -0.5), 0.3, 0.12, fc=CRIMSON, ec="none"))
c.text(0.6, -0.44, "first copy", ha="left", va="center", fontsize=8.8, color=TEXT2)
c.add_patch(mp.Rectangle((2.2, -0.5), 0.3, 0.12, fc=GREYF, ec="none"))
c.text(2.6, -0.44, "duplicate (same start, same end)", ha="left", va="center", fontsize=8.8,
       color=TEXT2)

# (d) blacklist --------------------------------------------------------------------------------------
d = fig.add_subplot(gs[1, 1]); d.set_xlim(0, 10); d.set_ylim(-0.6, 4.2); d.axis("off")
panel(d, "(d) ENCODE blacklist", "regions that light up in every experiment, even the input")
x = np.linspace(0, 10, 800)
def bump(c0, w, h): return h * np.exp(-0.5 * ((x - c0) / w) ** 2)
BL = (6.6, 8.0)
chip = 0.1 + bump(2.6, 0.18, 1.0) + bump(7.3, 0.35, 1.25) + 0.04 * rng.random(x.size)
inp = 0.1 + bump(7.3, 0.35, 1.15) + 0.04 * rng.random(x.size)
for y0, sig, lab, col in ((2.2, chip, "ChIP", tint(0.5)), (0.55, inp, "input", GREYF)):
    d.fill_between(x, y0, y0 + sig, color=col, lw=0)
    d.plot([0, 10], [y0, y0], color=LINE, lw=0.8)
    d.text(0.05, y0 + 0.95, lab, ha="left", va="center", fontsize=9.8, fontweight="semibold", color=TEXT)
d.add_patch(mp.Rectangle((BL[0], 0.45), BL[1] - BL[0], 3.3, fc="none", ec=DARK, lw=1.2, ls=(0, (3, 2))))
d.text(sum(BL) / 2, 3.85, "blacklist region", ha="center", va="bottom", fontsize=9.3, fontweight="semibold",
       color=DARK)
d.text(3.0, 3.0, "real peak\n(ChIP only)", ha="left", va="center", fontsize=9, color=CRIMSON,
       fontweight="semibold", linespacing=1.25)
d.text(5.0, -0.1, "repeats (satellite, centromere, rDNA) pile onto the few copies\n"
       "in the reference → remove before peak calling",
       ha="center", va="top", fontsize=9, color=TEXT2, linespacing=1.45)

header(a, "ChIP-seq QC", "Four checks before trusting a single peak", dy=40, x1=2.03)
fig.savefig(OUT + "chip_qc.png"); print("ok g24 chip qc")
