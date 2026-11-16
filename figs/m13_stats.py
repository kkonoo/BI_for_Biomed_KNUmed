from gallery import *
from m_sim import bh
from scipy import stats
rng = np.random.default_rng(13)

# Statistics for proteomics: (a) with n = 3 per group a protein's variance is a noisy guess, so limma shrinks it
# toward the typical value; (b) what that buys (simulated). Median normalisation is drawn in p3_norm.py.
fig, axes = plt.subplots(1, 2, figsize=(7.0, 4.2), gridspec_kw={"wspace": 0.3, "width_ratios": [1, 0.8]})
DARK = "#5E6670"; GREYF = "#B8BEC7"

# random draws kept from the former normalisation panel, so the simulation below is unchanged
base = rng.normal(23, 2, 1500); rng.normal(0, 0.35, (1500, 8))

# (b)/(c) 3 vs 3, 3,000 proteins, 10% truly change
P, n = 3000, 3
s_true = np.clip(np.sqrt(0.09 * 4 / rng.chisquare(4, P)), 0.1, 1.5)
eff = np.where(rng.random(P) < 0.1, rng.choice([-1, 1], P) * rng.uniform(0.6, 1.2, P), 0)
A = rng.normal(0, 1, (P, n)) * s_true[:, None]; B = rng.normal(0, 1, (P, n)) * s_true[:, None] + eff[:, None]
dif = B.mean(1) - A.mean(1)
s2 = (A.var(1, ddof=1) + B.var(1, ddof=1)) / 2; df = 2 * n - 2
d0, s0 = 4, np.median(s2)
s2m = (d0 * s0 + df * s2) / (d0 + df)
b = axes[0]
idx = rng.choice(P, 300, replace=False)
b.plot([0, 1.2], [0, 1.2], color=TEXT2, lw=0.8, ls=(0, (3, 2)))
b.scatter(s_true[idx], np.sqrt(s2[idx]), s=8, color=GREYF, lw=0, label="each protein alone")
b.scatter(s_true[idx], np.sqrt(s2m[idx]), s=8, color=CRIMSON, lw=0, label="moderated (limma)")
b.set_xlim(0, 1.0); b.set_ylim(0, 1.2)
b.annotate("SDs that look tiny\nby chance are pulled up", (0.2, 0.2), (0.45, 0.08), fontsize=8.3, color=CRIMSON,
           ha="left", va="center", arrowprops=dict(arrowstyle="-", color=CRIMSON, lw=0.7))
b.set_xlabel("true SD", color=TEXT2, fontsize=9); b.set_ylabel("estimated SD", color=TEXT2, fontsize=9)
for sp in ("left", "bottom"): b.spines[sp].set_color(LINE)
b.tick_params(colors=TEXT2, labelsize=8.3)
b.legend(frameon=False, fontsize=8.4, loc="upper left", labelcolor=TEXT, handletextpad=0.2, bbox_to_anchor=(-0.03, 1))
b.set_title("(a) Borrow variance", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
b.text(0, 1.03, "3 vs 3 samples, simulated", transform=b.transAxes, fontsize=9, color=TEXT2)

c = axes[1]
t = dif / np.sqrt(s2 * 2 / n); tm = dif / np.sqrt(s2m * 2 / n)
h1 = bh(2 * stats.t.sf(np.abs(t), df)) < 0.05
h2 = bh(2 * stats.t.sf(np.abs(tm), df + d0)) < 0.05
ntrue = (eff != 0).sum()
vals = [(h1 & (eff != 0)).sum(), (h2 & (eff != 0)).sum()]
c.bar([0, 1], vals, color=[GREYF, CRIMSON], width=0.55)
for x, v, h in zip([0, 1], vals, (h1, h2)):
    c.text(x, v + 3, f"{v}", ha="center", va="bottom", fontsize=11, fontweight="bold", color=TEXT)
    c.text(x, v + 3 + max(vals) * 0.09, f"+ {(h & (eff == 0)).sum()} false", ha="center", va="bottom", fontsize=8,
           color=TEXT2)
c.set_xticks([0, 1]); c.set_xticklabels(["t-test", "limma"], fontsize=9.5, color=TEXT)
c.tick_params(axis="x", length=0, pad=2); c.set_yticks([])
for sp in ("left", "top", "right"): c.spines[sp].set_visible(False)
c.spines["bottom"].set_color(LINE)
c.set_ylim(0, max(vals) * 1.4)
c.set_title("(b) Changes found", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
c.text(0, 1.03, f"of {ntrue}, at 5% FDR", transform=c.transAxes, fontsize=9, color=TEXT2)
header(b, "Proteomics · statistics", "Same logic as RNA-seq: borrow strength across proteins",
       dy=40, x1=2.05)
import os
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "assets", "img", "w6") + "/"
fig.savefig(OUT + "prot_stats.png"); print("ok m13", vals, h1.sum(), h2.sum(), ntrue)
