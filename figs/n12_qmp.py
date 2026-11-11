from gallery import *

# Quantitative microbiome profiling: share × total load = amount. Example (illustrative numbers): after an
# antibiotic the total load falls ten-fold; one taxon that survives better rises in share but falls in amount.
fig, axes = plt.subplots(1, 3, figsize=(9.8, 3.9), gridspec_kw={"wspace": 0.5})
DARK = "#5E6670"; GREYF = "#B8BEC7"
share = np.array([0.10, 0.60]); load = np.array([10.0, 1.0]); amount = share * load
panels = [("(a) Share", "what sequencing reports", share, "share of reads", lambda v: f"{v:.0%}", CRIMSON),
          ("(b) Total load", "flow cytometry or qPCR", load, "10¹⁰ cells per gram", lambda v: f"{v:g}", DARK),
          ("(c) Share × load", "the amount of this taxon", amount, "10¹⁰ cells per gram", lambda v: f"{v:g}", DARK)]
ops = ["×", "="]
for k, (ax, (t, sub, v, yl, fmt, col)) in enumerate(zip(axes, panels)):
    ax.bar([0, 1], v, width=0.55, color=[GREYF, col])
    for x, val in enumerate(v):
        ax.text(x, val + v.max() * 0.03, fmt(val), ha="center", va="bottom", fontsize=9.5, fontweight="semibold",
                color=TEXT)
    r = v[1] / v[0]
    ax.text(0.5, v.max() * 1.22, f"×{r:g}", ha="center", va="center", fontsize=11, fontweight="bold",
            color=CRIMSON if r > 1 else DARK)
    ax.set_ylim(0, v.max() * 1.35)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["before", "after\nantibiotic"], fontsize=8.8, color=TEXT)
    ax.set_yticks([]); ax.set_ylabel(yl, color=TEXT2, fontsize=8.6)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    ax.spines["left"].set_color(LINE); ax.spines["bottom"].set_color(LINE)
    ax.tick_params(axis="x", length=0)
    ax.set_title(t, loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
    ax.text(0, 1.03, sub, transform=ax.transAxes, fontsize=9, color=TEXT2)
    if k < 2:
        ax.text(1.27, 0.45, ops[k], transform=ax.transAxes, ha="center", va="center", fontsize=18, color=TEXT2)
fig.text(0.125, -0.07, "Illustrative numbers. The taxon 'rises 6-fold' in share while its cells fall by 40%: "
         "only a measured total turns shares into amounts.", fontsize=8.8, color=TEXT2)
header(axes[0], "Metagenomics · quantitative profiling", "A rising share can hide a falling amount", dy=40, x1=3.75)
fig.savefig(OUT + "qmp.png"); print("ok n12", amount)
