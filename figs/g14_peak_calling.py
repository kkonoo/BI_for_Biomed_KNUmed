from gallery import *
from math import lgamma, log, exp

rng = np.random.default_rng(7)

def bump(x, c, w, h):
    return h * np.exp(-0.5 * ((x - c) / w) ** 2)

def pois_sf(k, lam):
    """P(X >= k) for X ~ Poisson(lam), summed in log space."""
    terms = [exp(i * log(lam) - lam - lgamma(i + 1)) for i in range(k, k + 400)]
    return sum(terms)

fig = plt.figure(figsize=(11.5, 9.6))
gs = fig.add_gridspec(2, 1, height_ratios=[1.15, 1.0], hspace=0.36)
top = gs[0].subgridspec(3, 2, width_ratios=[2.55, 1], height_ratios=[1, 0.8, 0.32],
                        hspace=0.18, wspace=0.16)
bot = gs[1].subgridspec(3, 2, height_ratios=[1, 0.36, 0.36], hspace=0.12, wspace=0.14)

# ① one window vs its local background ------------------------------------------------------------
x = np.linspace(0, 20, 801)                       # kb
inp = 2.2 + 0.25 * np.sin(x / 2.3) + bump(x, 14.2, 0.45, 3.4) + rng.normal(0, 0.28, x.size).clip(-0.6, 0.6)
chip = 2.4 + bump(x, 6.0, 0.32, 15.5) + bump(x, 14.2, 0.45, 4.6) + rng.normal(0, 0.45, x.size).clip(-0.9, 0.9)
inp = np.clip(inp, 0.2, None); chip = np.clip(chip, 0.2, None)

a1 = fig.add_subplot(top[0, 0]); a2 = fig.add_subplot(top[1, 0], sharex=a1); a3 = fig.add_subplot(top[2, 0], sharex=a1)
win = (5.6, 6.4)
for a in (a1, a2):
    a.axvspan(*win, color=tint(0.10), zorder=0, lw=0)
a1.fill_between(x, chip, color=tint(0.32), lw=0)
a1.fill_between(x, chip, where=(x > win[0]) & (x < win[1]), color=CRIMSON, lw=0)
a2.fill_between(x, inp, color="#B8BEC7", lw=0)
for a, lab in ((a1, "ChIP"), (a2, "input")):
    a.set_yticks([]); a.spines["left"].set_visible(False); a.spines["bottom"].set_color(LINE)
    a.tick_params(axis="x", colors=TEXT2, labelbottom=False, length=0)
    a.text(-0.3, a.get_ylim()[1] * 0.5 if a is a2 else 8, lab, ha="right", va="center", fontsize=11,
           fontweight="semibold", color=TEXT)
a1.set_ylim(0, 19); a2.set_ylim(0, 8.5)
a1.text(6.55, 16.4, "window: 48 reads", ha="left", va="center", fontsize=10, color=CRIMSON, fontweight="semibold")
a1.text(14.2, 8.0, "also high in input", ha="center", va="bottom", fontsize=9.5, color=TEXT2)

# local lambda: background measured in the input around the window
for y, span, lab in ((7.6, 1, "1 kb"), (6.7, 5, "5 kb"), (5.8, 10, "10 kb")):
    a2.plot([6 - span / 2, 6 + span / 2], [y, y], color=TEXT2, lw=1.1, solid_capstyle="butt")
    for e in (6 - span / 2, 6 + span / 2):
        a2.plot([e, e], [y - 0.22, y + 0.22], color=TEXT2, lw=1.1)
    a2.text(6 + span / 2 + 0.15, y, lab, ha="left", va="center", fontsize=8.8, color=TEXT2)
a3.text(0, -0.75, r"background from the input:  $\lambda_{local}$ = max($\lambda_{genome}$, $\lambda_{1\,kb}$, $\lambda_{5\,kb}$, $\lambda_{10\,kb}$)  = 6 reads expected in the window",
        ha="left", va="center", fontsize=9.8, color=TEXT)

# called peaks lane
a3.set_ylim(0, 1); a3.axis("off")
a3.add_patch(mp.Rectangle((win[0], 0.25), win[1] - win[0], 0.5, fc=CRIMSON, ec="none"))
a3.text(win[1] + 0.2, 0.5, "peak called", ha="left", va="center", fontsize=9.5, color=CRIMSON, fontweight="semibold")
a3.add_patch(mp.Rectangle((13.5, 0.25), 1.4, 0.5, fc="white", ec=LINE, lw=1, ls=(0, (3, 2))))
a3.text(15.1, 0.5, "not called: ChIP ≈ input", ha="left", va="center", fontsize=9.5, color=TEXT2)
a3.text(-0.3, 0.5, "calls", ha="right", va="center", fontsize=11, fontweight="semibold", color=TEXT)
a1.set_xlim(0, 20)
a2.tick_params(axis="x", labelbottom=False)
a2.text(20, -1.2, "20 kb (schematic)", ha="right", va="top", fontsize=8.8, color=TEXT2)

# Poisson inset
ai = fig.add_subplot(top[:, 1])
lam, k = 6, 48
ks = np.arange(0, 56)
pmf = np.array([exp(i * log(lam) - lam - lgamma(i + 1)) for i in ks])
ai.bar(ks, pmf, width=0.85, color="#B8BEC7", lw=0)
ai.axvline(k, color=CRIMSON, lw=2)
p = pois_sf(k, lam)
ai.text(k - 1.2, pmf.max() * 0.92, f"observed = {k}", ha="right", va="center", fontsize=10,
        color=CRIMSON, fontweight="semibold")
ai.text(k - 1.2, pmf.max() * 0.80, f"P(X ≥ {k}) ≈ " + r"$10^{%d}$" % round(np.log10(p)), ha="right", va="center",
        fontsize=10, color=CRIMSON)
ai.text(13, pmf.max() * 0.58, "Poisson(λ = 6):\nreads the background\npredicts", ha="left", va="center",
        fontsize=9.3, color=TEXT2, linespacing=1.35)
ai.text(11, pmf.max() * 0.25, "→ q value\n   (all windows tested)\n→ q < cutoff = peak",
        ha="left", va="center", fontsize=9.3, color=TEXT, linespacing=1.4)
ai.set_xlim(-1, 55); ai.set_yticks([]); ai.spines["left"].set_visible(False)
ai.spines["bottom"].set_color(LINE); ai.tick_params(colors=TEXT2)
ai.set_xlabel("reads in the window", color=TEXT2, fontsize=9.5)
a1.set_title("① Is this window enriched over the local background?", loc="left", fontsize=12.5,
             fontweight="semibold", color=TEXT, pad=8)
header(a1, "ChIP-seq · peak calling (MACS2)", "A peak = more reads than the input background predicts", dy=24, x1=1.42)

# ② narrow vs broad ----------------------------------------------------------------------------------
xs = np.linspace(0, 30, 1201)
sharp = 1.0 + bump(xs, 6, 0.28, 12) + bump(xs, 17, 0.28, 9.5) + bump(xs, 24, 0.3, 11) + rng.normal(0, 0.35, xs.size).clip(-0.7, 0.7)
xb = np.linspace(0, 120, 1201)
dom = (xb > 25) & (xb < 95)
broad = 1.0 + 2.3 * dom * (1 + 0.35 * np.sin(xb / 3.1)) + rng.normal(0, 0.55, xb.size).clip(-1, 1)
broad = np.clip(broad, 0.2, None); sharp = np.clip(sharp, 0.2, None)

def track(ax, xx, yy, xmax, label):
    ax.fill_between(xx, yy, color=tint(0.45), lw=0)
    ax.set_ylim(0, 14.5); ax.set_xlim(0, xmax); ax.set_yticks([]); ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_color(LINE); ax.tick_params(axis="x", labelbottom=False, length=0)
    ax.text(0, ax.get_ylim()[1] * 1.0, label, ha="left", va="top", fontsize=10, color=TEXT2)

def lane(ax, xmax, segs, fc, ec, label, lc):
    ax.set_xlim(0, xmax); ax.set_ylim(0, 1); ax.axis("off")
    for a, b in segs:
        ax.add_patch(mp.Rectangle((a, 0.05), b - a, 0.42, fc=fc, ec=ec, lw=0.9))
    ax.text(0, 0.8, label, ha="left", va="center", fontsize=9.3, color=lc, fontweight="semibold")

b1 = fig.add_subplot(bot[0, 0]); track(b1, xs, sharp, 30, "H3K4me3-like: sharp peaks at TSSs")
b1.text(30, 13.8, "~30 kb (schematic)", ha="right", va="top", fontsize=8.8, color=TEXT2)
b1.set_title("② Narrow vs broad mode", loc="left", fontsize=12.5, fontweight="semibold", color=TEXT, pad=8)
lane(fig.add_subplot(bot[1, 0]), 30, [(5.3, 6.7), (16.3, 17.7), (23.3, 24.7)], CRIMSON, "none",
     "narrow mode: one call per peak", CRIMSON)
ax_note = fig.add_subplot(bot[2, 0]); ax_note.axis("off")
ax_note.text(0, 0.5, "TFs, H3K4me3, H3K27ac — peaks a few hundred bp wide", ha="left", va="center",
             fontsize=9.3, color=TEXT2)

b2 = fig.add_subplot(bot[0, 1]); track(b2, xb, broad, 120, "H3K27me3-like: low, wide enrichment")
frag = [(28, 31), (37, 39.5), (47, 50), (58, 60), (66, 69.5), (79, 81), (88, 91)]
lane(fig.add_subplot(bot[1, 1]), 120, frag, "white", TEXT2, "narrow mode: fragments", TEXT2)
lane(fig.add_subplot(bot[2, 1]), 120, [(25, 95)], CRIMSON, "none", "broad mode: one domain", CRIMSON)
b2.text(120, 13.8, "~100 kb (schematic)", ha="right", va="top", fontsize=8.8, color=TEXT2)

fig.savefig(OUT + "peak_calling.png"); print("ok g14 peak calling", p)
