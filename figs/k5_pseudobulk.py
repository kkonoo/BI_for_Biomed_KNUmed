from gallery import *
from scipy import stats
rng = np.random.default_rng(5)

# Per-cell tests vs pseudobulk under the null (no condition effect, real donor-to-donor variation):
# the per-cell false-positive rate climbs with the number of cells; pseudobulk stays at 5%
fig = plt.figure(figsize=(9.8, 4.7))
gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.05], wspace=0.25, left=0.02, right=0.98, top=0.8, bottom=0.14)
DARK = "#5E6670"; GREYF = "#B8BEC7"

# (a) the data structure
a = fig.add_subplot(gs[0]); a.set_xlim(0, 10); a.set_ylim(0, 10); a.axis("off")
a.set_title("(a) Who is the replicate?", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
a.text(0, 1.03, "3 patients vs 3 controls, 900 cells each", transform=a.transAxes, fontsize=9, color=TEXT2)
dm = [0.5, -0.2, 0.25, -0.35, 0.1, -0.15]
for k in range(6):
    x = 0.9 + k * 1.55 + (0.35 if k >= 3 else 0)
    c = CRIMSON if k < 3 else DARK
    pts = rng.normal([x, 6.6 + dm[k]], [0.28, 0.55], (40, 2))
    a.scatter(pts[:, 0], pts[:, 1], s=6, color=tint(0.45) if k < 3 else GREYF, lw=0)
    a.add_patch(mp.Rectangle((x - 0.38, 2.6), 0.76, 0.55, fc=c, ec="none"))
    a.text(x, 2.875, f"{'P' if k < 3 else 'C'}{k % 3 + 1}", ha="center", va="center", fontsize=8.6,
           color="white", fontweight="bold")
    a.annotate("", (x, 3.25), (x, 4.6), arrowprops=dict(arrowstyle="-|>", color=LINE, lw=0.9, mutation_scale=7))
a.text(0.1, 9.25, "per-cell test: every dot counts as a sample", ha="left", fontsize=9, color=TEXT)
a.text(9.9, 9.25, "n = 5,400", ha="right", fontsize=9.5, color=CRIMSON, fontweight="bold")
a.text(0.1, 3.95, "sum counts per donor", ha="left", va="center", fontsize=8.5, color=TEXT2, zorder=5,
       bbox=dict(fc="white", ec="none", pad=2))
a.text(0.1, 1.75, "pseudobulk: one sample per donor → DESeq2 / edgeR", ha="left", fontsize=9, color=TEXT)
a.text(9.9, 1.75, "n = 6", ha="right", fontsize=9.5, color=TEXT, fontweight="bold")
a.text(0.1, 0.65, "cells from one person share genetics, drugs,\nbatch and handling → they are not independent",
       ha="left", va="center", fontsize=8.5, color=TEXT2, linespacing=1.25)

# (b) false-positive rate vs cells per donor, simulated with no true difference
b = fig.add_subplot(gs[1])
ns = [10, 30, 100, 300, 1000]; R = 400; fp_cell, fp_pb = [], []
for n in ns:
    hc = hp = 0
    for _ in range(R):
        mu = rng.normal(0, 0.25, 6)                       # donor effects; no patient-vs-control effect
        x = [rng.normal(m, 1.0, n) for m in mu]
        p1 = stats.mannwhitneyu(np.concatenate(x[:3]), np.concatenate(x[3:])).pvalue
        p2 = stats.ttest_ind([v.mean() for v in x[:3]], [v.mean() for v in x[3:]]).pvalue
        hc += p1 < 0.05; hp += p2 < 0.05
    fp_cell.append(hc / R); fp_pb.append(hp / R)
b.plot(ns, fp_cell, "o-", color=CRIMSON, lw=2, ms=5)
b.plot(ns, fp_pb, "o-", color=DARK, lw=2, ms=5)
b.axhline(0.05, color=LINE, lw=0.9, ls=(0, (4, 3)))
b.text(ns[-1], fp_cell[-1] + 0.05, f"per-cell Wilcoxon  {100 * fp_cell[-1]:.0f}%", ha="right", fontsize=9,
       color=CRIMSON, fontweight="semibold")
b.text(ns[-1], fp_pb[-1] + 0.05, f"pseudobulk  {100 * fp_pb[-1]:.0f}%  (= the nominal 5%)", ha="right", fontsize=9, color=DARK,
       fontweight="semibold")
b.set_xscale("log"); b.set_ylim(0, 1); b.set_xticks(ns); b.set_xticklabels([str(n) for n in ns])
b.set_yticks([0, 0.25, 0.5, 0.75, 1]); b.set_yticklabels(["0", "25%", "50%", "75%", "100%"])
b.set_xlabel("cells per donor", color=TEXT2, fontsize=9.3); b.set_ylabel("false-positive rate", color=TEXT2, fontsize=9.3)
for sp in ("left", "bottom"): b.spines[sp].set_color(LINE)
b.tick_params(colors=TEXT2, labelsize=8.3); b.minorticks_off()
b.set_title("(b) Simulated: no real difference", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
b.text(0, 1.03, f"share of {R} runs called significant at p < 0.05", transform=b.transAxes, fontsize=9, color=TEXT2)
header(a, "Single-cell RNA-seq · statistics", "More cells make a per-cell test more wrong, not more powerful",
       dy=40, x1=2.2)
fig.savefig(OUT + "pseudobulk.png"); print("ok k5 pseudobulk", fp_cell, fp_pb)
