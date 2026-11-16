from gallery import *
import os
rng = np.random.default_rng(13)
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "assets", "img", "w6") + "/"

# Normalisation of protein intensities (simulated): (a) raw intensities are extremely skewed;
# (b) on the log2 scale they are roughly symmetric and a fold change is a fixed distance;
# (c) samples loaded with different amounts shift as a whole; median centring lines them up.
fig, axes = plt.subplots(1, 3, figsize=(10.4, 4.1), gridspec_kw={"wspace": 0.34, "width_ratios": [1, 1, 1.25]})
DARK = "#5E6670"; GREYF = "#B8BEC7"

def style(a, title, sub):
    for sp in ("left", "bottom"): a.spines[sp].set_color(LINE)
    a.tick_params(colors=TEXT2, labelsize=8.3)
    a.set_title(title, loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
    a.text(0, 1.03, sub, transform=a.transAxes, fontsize=9, color=TEXT2)

L = rng.normal(23, 2, 4000)            # log2 intensities of 4,000 proteins
R = 2 ** L

a = axes[0]
a.hist(R / 1e7, bins=np.linspace(0, 120, 61), color=GREYF, lw=0)
a.set_xlim(0, 120); a.set_yticks([]); a.spines["left"].set_visible(False)
a.set_xlabel("intensity (× 10⁷)", color=TEXT2, fontsize=9)
a.annotate("most proteins are\nsquashed near 0", (3, 1700), (25, 2300), fontsize=8.4, color=TEXT2, ha="left",
           va="center", arrowprops=dict(arrowstyle="-", color=TEXT2, lw=0.7))
a.annotate("a few abundant ones\nstretch the axis", (100, 25), (55, 900), fontsize=8.4, color=TEXT2, ha="left",
           va="center", arrowprops=dict(arrowstyle="-", color=TEXT2, lw=0.7))
style(a, "(a) Raw intensities", "4,000 proteins, one sample")

b = axes[1]
b.hist(L, bins=np.linspace(15, 31, 49), color=DARK, lw=0)
b.set_xlim(15, 31); b.set_yticks([]); b.spines["left"].set_visible(False)
b.set_xlabel("log₂ intensity", color=TEXT2, fontsize=9)
cnt, _ = np.histogram(L, bins=np.linspace(15, 31, 49)); y0 = cnt.max() * 0.2
for x1 in (16.9, 28.1):
    b.plot([x1, x1 + 1], [y0, y0], color=CRIMSON, lw=1.6, solid_capstyle="butt")
    for xe in (x1, x1 + 1):
        b.plot([xe, xe], [y0 - cnt.max() * 0.025, y0 + cnt.max() * 0.025], color=CRIMSON, lw=1.6)
    b.text(x1 + 0.5, y0 + cnt.max() * 0.04, "×2", fontsize=8.6, color=CRIMSON, ha="center", va="bottom")
b.text(30.8, cnt.max() * 0.95, "a doubling is +1\nanywhere on\nthe axis", fontsize=8.4, color=CRIMSON,
       ha="right", va="top")
style(b, "(b) After log₂", "same proteins")

c = axes[2]
off = np.array([0.0, 0.8, -0.6, 0.3, 1.1, -0.9, 0.5, -0.3])
D = L[:1500, None] + off[None, :] + rng.normal(0, 0.35, (1500, 8))
Dn = D - np.median(D, 0) + np.median(D)
for M, x0, col in ((D, 0, GREYF), (Dn, 9, DARK)):
    bp = c.boxplot(M, positions=np.arange(8) + x0, widths=0.6, patch_artist=True, showfliers=False,
                   medianprops=dict(color=CRIMSON, lw=1.6), whiskerprops=dict(color=LINE), capprops=dict(color=LINE),
                   boxprops=dict(lw=0))
    for p in bp["boxes"]: p.set_facecolor(col)
c.set_xticks([3.5, 12.5]); c.set_xticklabels(["before", "after"], fontsize=9.5, color=TEXT)
c.tick_params(axis="x", length=0)
c.set_ylabel("log₂ intensity", color=TEXT2, fontsize=9)
style(c, "(c) Median normalisation", "8 samples; red line = median")
c.tick_params(axis="x", length=0)

header(a, "Proteomics · normalisation", "Log-transform, then line the samples up", dy=40, x1=3.55)
fig.savefig(OUT + "prot_norm.png"); print("ok p3 norm")
