from gallery import *
rng = np.random.default_rng(15)

# The same DE result shown two ways: volcano (significance vs effect) and MA (effect vs expression)
G = 12000
mu = np.exp(rng.normal(3.6, 2.0, G))
lfc = rng.normal(0, 0.32, G)
idx = rng.choice(G, 420, replace=False)
lfc[idx] += rng.choice([-1, 1], 420) * rng.gamma(2.4, 0.55, 420)
se = 0.1 + 1.8 / np.sqrt(mu + 8)
z = lfc / se
from scipy.stats import norm
p = np.clip(2 * norm.sf(np.abs(z)), 1e-300, 1)
o_ = np.argsort(p); bh = p[o_] * G / np.arange(1, G + 1)          # Benjamini–Hochberg
bh = np.minimum.accumulate(bh[::-1])[::-1]
padj = np.empty(G); padj[o_] = np.clip(bh, 1e-300, 1)
sig = (padj < 0.05) & (np.abs(lfc) > 1)
DOWN = "#5E6670"; NS = "#CDD1D7"
up, dn = sig & (lfc > 0), sig & (lfc < 0)
col = np.where(up, CRIMSON, np.where(dn, DOWN, NS))
order = np.argsort(sig)                         # draw significant genes on top

fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.8), gridspec_kw={"wspace": 0.22})
def style(ax, title, sub):
    for sp in ("left", "bottom"): ax.spines[sp].set_color(LINE)
    ax.tick_params(colors=TEXT2, labelsize=8.5)
    ax.set_title(title, loc="left", fontsize=12.5, fontweight="semibold", color=TEXT, pad=24)
    ax.text(0, 1.03, sub, transform=ax.transAxes, ha="left", va="bottom", fontsize=9.3, color=TEXT2)

# (a) volcano
ax = axes[0]
yv = -np.log10(padj); cap = yv > 40; yv = np.minimum(yv, 40)
o1 = order[~cap[order]]
ax.scatter(lfc[o1], yv[o1], s=5, c=col[o1], alpha=0.75, lw=0)
ax.scatter(lfc[cap], yv[cap], s=14, c=col[cap], marker="^", lw=0)        # beyond the axis
ax.axhline(-np.log10(0.05), color=TEXT2, ls=(0, (3, 2)), lw=1)
for v in (-1, 1): ax.axvline(v, color=TEXT2, ls=(0, (1, 2)), lw=1)
ax.set_xlim(-6, 6); ax.set_ylim(0, 46)
ax.set_xlabel("log₂ fold change", color=TEXT2, fontsize=9)
ax.set_ylabel("−log₁₀(adjusted p)", color=TEXT2, fontsize=9)
yt = ax.get_ylim()[1]
ax.text(-5.8, yt * 0.97, f"{dn.sum()} down", ha="left", va="top", fontsize=10, color=DOWN, fontweight="bold")
ax.text(5.8, yt * 0.97, f"{up.sum()} up", ha="right", va="top", fontsize=10, color=CRIMSON, fontweight="bold")
ax.text(5.8, -np.log10(0.05) + 0.4, "padj = 0.05", ha="right", va="bottom", fontsize=8.5, color=TEXT2)
ax.text(0, yt * 0.97, "|log₂FC| = 1", ha="center", va="top", fontsize=8.5, color=TEXT2,
        bbox=dict(fc="white", ec="none", pad=1))
ax.text(0, yt * 0.72, "⚠ both lines are\nconventions, not facts", ha="center", va="center", fontsize=8.8,
        color=TEXT, linespacing=1.35, bbox=dict(fc="white", ec="none", pad=1.5))
style(ax, "(a) Volcano plot", "significance vs effect size · ▲ = off the scale")

# (b) MA
ax = axes[1]
ax.scatter(mu[order], lfc[order], s=5, c=col[order], alpha=0.75, lw=0)
ax.set_xscale("log"); ax.axhline(0, color=TEXT2, lw=1)
ax.set_ylim(-6, 6); ax.set_xlim(0.3, 2e5)
ax.set_xlabel("mean normalised count (expression)", color=TEXT2, fontsize=9)
ax.set_ylabel("log₂ fold change", color=TEXT2, fontsize=9)
ax.annotate("low counts: huge fold changes,\nbut too noisy to be significant", (1.2, -3.0), (40, -5.0),
            fontsize=8.8, color=TEXT2, ha="left", va="center", linespacing=1.35,
            arrowprops=dict(arrowstyle="-|>", color=TEXT2, lw=1, connectionstyle="arc3,rad=-0.2"))
ax.text(ax.get_xlim()[1] * 0.8, 5.6, "hits should not pile up\nat low counts", ha="right", va="top",
        fontsize=8.8, color=TEXT, linespacing=1.35)
style(ax, "(b) MA plot", "effect size vs expression · the better diagnostic")
for k, (lab, c) in enumerate((("up", CRIMSON), ("down", DOWN), ("not significant", "#A9AFB8"))):
    ax.text((0.03, 0.15, 0.3)[k], 0.96, "● " + lab, transform=ax.transAxes, ha="left", va="top", fontsize=8.8,
            color=c, fontweight="semibold")

header(axes[0], "Differential expression · reading the result",
       "Same genes, two views: volcano for significance, MA for diagnosis", dy=44, x1=2.22)
fig.savefig(OUT + "volcano_ma.png"); print("ok h4")
