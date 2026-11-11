from gallery import *
from scipy import stats
rng = np.random.default_rng(14)

# Confounding by medication (simulated). Escherichia abundance depends on metformin, not on diabetes.
# (a) cases vs controls looks like a diabetes signature; (b) split the cases by treatment and it belongs to the drug.
n_c, n_d = 250, 250
met = rng.random(n_d) < 0.6
ctrl = rng.normal(0, 1, n_c)
t2d = rng.normal(0, 1, n_d) + 1.3 * met
fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.2), gridspec_kw={"wspace": 0.25, "width_ratios": [1, 1.4]})
DARK = "#5E6670"; GREYF = "#B8BEC7"

def strip(ax, x, v, col, lab):
    ax.scatter(x + rng.uniform(-0.2, 0.2, len(v)), v, s=7, color=col, lw=0, alpha=0.65)
    m = np.median(v)
    ax.plot([x - 0.3, x + 0.3], [m, m], color=TEXT, lw=2)
    ax.text(x, -3.6, lab, ha="center", va="top", fontsize=9, color=TEXT, linespacing=1.1)

def style(ax, t, sub):
    ax.set_title(t, loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
    ax.text(0, 1.03, sub, transform=ax.transAxes, fontsize=9, color=TEXT2)
    ax.set_xticks([]); ax.set_ylim(-3.4, 4.8)
    for sp in ("top", "right", "bottom"): ax.spines[sp].set_visible(False)
    ax.spines["left"].set_color(LINE); ax.tick_params(colors=TEXT2, labelsize=8.3)

a = axes[0]
strip(a, 0, ctrl, GREYF, f"controls\n(n = {n_c})")
strip(a, 1, t2d, CRIMSON, f"diabetes\n(n = {n_d})")
p = stats.mannwhitneyu(t2d, ctrl).pvalue
sup = str(int(np.ceil(np.log10(p)))).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
a.text(0.5, 4.3, f"p < 10{sup}", ha="center", va="center", fontsize=10, fontweight="bold", color=CRIMSON)
a.set_xlim(-0.6, 1.6); a.set_ylabel("Escherichia (log abundance)", color=TEXT2, fontsize=9)
style(a, "(a) Cases vs controls", "looks like a diabetes signature")

b = axes[1]
strip(b, 0, ctrl, GREYF, "controls")
strip(b, 1, t2d[~met], GREYF, f"diabetes,\nno metformin")
strip(b, 2, t2d[met], CRIMSON, f"diabetes,\non metformin")
p1 = stats.mannwhitneyu(t2d[~met], ctrl).pvalue
b.text(0.5, 4.3, f"p = {p1:.2f}", ha="center", va="center", fontsize=10, fontweight="bold", color=DARK)
b.annotate("", (2, 3.7), (1, 3.7), arrowprops=dict(arrowstyle="-|>", color=CRIMSON, lw=1.2, mutation_scale=9))
b.text(1.5, 4.3, "the difference is here", ha="center", va="center", fontsize=9.5, fontweight="bold", color=CRIMSON)
b.set_xlim(-0.6, 2.6)
style(b, "(b) Split the cases by treatment", f"{met.mean():.0%} of the cases take metformin; no control does")
header(a, "Metagenomics · confounding", "A 'diabetes signature' that belonged to the drug", dy=40, x1=2.5)
fig.text(0.125, -0.08, "Simulated: Escherichia responds to metformin only. The real case is Forslund et al., Nature 2015.",
         fontsize=8.8, color=TEXT2)
fig.savefig(OUT + "metformin.png"); print("ok n14 p", p, "p untreated", round(p1, 3), "met", met.mean())
