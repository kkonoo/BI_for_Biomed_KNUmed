from gallery import *
from scipy.stats import norm

# Power depends on replicates (left); extra depth saturates, extra replicates do not (right)
fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.6), gridspec_kw={"wspace": 0.22})
DARK = "#5E6670"; GREYF = "#B8BEC7"

def style(ax, title, sub):
    for sp in ("left", "bottom"): ax.spines[sp].set_color(LINE)
    ax.tick_params(colors=TEXT2, labelsize=8.5)
    ax.set_title(title, loc="left", fontsize=12.5, fontweight="semibold", color=TEXT, pad=24)
    ax.text(0, 1.03, sub, transform=ax.transAxes, ha="left", va="bottom", fontsize=9.3, color=TEXT2)

# (a) power vs n
ax = axes[0]
n = np.arange(2, 31)
for lfc, c, lab, ly in [(0.5, GREYF, "1.4×  (log₂FC = 0.5)", None), (1.0, tint(0.5), "2×  (log₂FC = 1)", None),
                        (2.0, CRIMSON, "4×  (log₂FC = 2)", None)]:
    se = np.sqrt(2) * 0.4 / np.sqrt(n)
    zcrit = norm.isf(0.05 / 20000 / 2)
    pw = norm.sf(zcrit - (lfc * np.log(2) * 0.7) / se)
    ax.plot(n, pw, color=c, lw=2.2)
    xs, ys, va = {0.5: (19, 0.035, "bottom"), 1.0: (23, 0.36, "bottom"), 2.0: (16, 0.92, "top")}[lfc]
    ax.text(xs, ys, lab, ha="right" if lfc == 1.0 else "left", va=va, fontsize=8.8, color=TEXT2 if c == GREYF else (CRIMSON if lfc == 2 else tint(0.8)),
            fontweight="semibold")
ax.axhline(0.8, color=LINE, ls=(0, (3, 2)), lw=1)
ax.text(30, 0.83, "80% power", ha="right", va="bottom", fontsize=8.8, color=TEXT2)
ax.axvline(3, color=DARK, ls=(0, (3, 2)), lw=1.2)
ax.text(3.5, 0.6, "n = 3", ha="left", va="center", fontsize=9.5, color=DARK, fontweight="semibold")
ax.set_xlabel("biological replicates per group", color=TEXT2, fontsize=9)
ax.set_ylabel("power to detect a gene", color=TEXT2, fontsize=9)
ax.set_ylim(0, 1.02); ax.set_xlim(2, 30)
style(ax, "(a) Power depends on n", "lines = effect size of the gene")

# (b) depth vs replicates
ax = axes[1]
depth = np.array([2, 5, 10, 20, 40, 80])
base = np.array([0.30, 0.52, 0.68, 0.76, 0.79, 0.80])
ax.axvspan(20, 80, color=FAINT, lw=0, zorder=0)
for r, c in zip([2, 3, 5, 10], [GREYF, DARK, tint(0.5), CRIMSON]):
    scale = {2: 0.35, 3: 0.55, 5: 0.80, 10: 1.0}[r]
    y = base * scale
    ax.plot(depth, y, marker="o", ms=4.5, color=c, lw=2)
    ax.text(92, y[-1], f"n = {r}", ha="left", va="center", fontsize=9, fontweight="semibold",
            color=TEXT2 if c == GREYF else (CRIMSON if c != DARK else DARK))
ax.set_xscale("log"); ax.set_xticks(depth); ax.set_xticklabels([str(d) for d in depth])
ax.minorticks_off()
ax.set_xlim(1.7, 140)
ax.set_xlabel("sequencing depth (million reads per sample)", color=TEXT2, fontsize=9)
ax.set_ylabel("relative sensitivity", color=TEXT2, fontsize=9)
ax.set_ylim(0, 0.9)
ax.text(40, 0.06, "20–30 M reads is enough;\nmore depth buys almost nothing", ha="center", va="bottom", fontsize=8.8,
        color=TEXT2, linespacing=1.4)
ax.annotate("", (2.15, 0.84), (2.15, 0.36), arrowprops=dict(arrowstyle="-|>", color=CRIMSON, lw=1.3))
ax.text(2.45, 0.86, "more replicates\nlift every curve", ha="left", va="top", fontsize=8.8, color=CRIMSON,
        fontweight="semibold", linespacing=1.3)
style(ax, "(b) Depth saturates, replicates do not", "spend the money on more samples")

header(axes[0], "Differential expression · study size",
       "Power comes from replicates, not from deeper sequencing", dy=44, x1=2.22)
fig.savefig(OUT + "power_replicates.png"); print("ok h5")
