from gallery import *
rng = np.random.default_rng(5)

# Compositional data. 30 taxa; in cases one taxon blooms (5% → 50% of the community, a 19-fold rise in cells);
# the other 29 do not change at all. (a) cells; (b) proportions; (c) how one unchanged taxon looks by method.
ntax = 30
others = np.exp(rng.normal(0, 0.8, ntax - 1)); others = 95 * others / others.sum()   # 95 units of cells
ctrl = np.r_[5.0, others]; case = np.r_[95.0, others]                               # bloomer: 5 → 95 units
fig, axes = plt.subplots(1, 3, figsize=(9.8, 4.3), gridspec_kw={"wspace": 0.55, "width_ratios": [1, 1, 1.25]})
DARK = "#5E6670"; GREYF = "#B8BEC7"

def style(ax, t, sub):
    ax.set_title(t, loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
    ax.text(0, 1.03, sub, transform=ax.transAxes, fontsize=9, color=TEXT2)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"): ax.spines[sp].set_color(LINE)
    ax.tick_params(colors=TEXT2, labelsize=8.3)

def stacked(ax, cols):
    for x, v in enumerate(cols):
        ax.bar(x, v[0], color=CRIMSON, width=0.55, edgecolor="white", lw=0.6)
        bot = v[0]
        for k, u in enumerate(v[1:]):
            ax.bar(x, u, bottom=bot, color=GREYF if k % 2 else "#CDD2D9", width=0.55, edgecolor="white", lw=0.3)
            bot += u
    ax.set_xticks([0, 1]); ax.set_xticklabels(["control", "case"], fontsize=9.2, color=TEXT)
    ax.tick_params(axis="x", length=0)

a = axes[0]
stacked(a, [ctrl, case]); a.set_ylabel("bacterial cells (arbitrary units)", color=TEXT2, fontsize=9)
a.text(1.33, 95 + 95 / 2, "the other 29:\nunchanged", ha="left", va="center", fontsize=8.5, color=TEXT2, linespacing=1.15)
a.text(1.33, 95 / 2, "one taxon\n×19", ha="left", va="center", fontsize=8.5, color=CRIMSON, fontweight="semibold",
       linespacing=1.15)
style(a, "(a) Cells", "what happened in the gut")

b = axes[1]
stacked(b, [ctrl / ctrl.sum(), case / case.sum()])
b.set_yticks([0, 0.5, 1]); b.set_yticklabels(["0", "50%", "100%"])
b.text(1.33, 0.75, "the other 29:\n95% → 50%", ha="left", va="center", fontsize=8.5, color=TEXT2, linespacing=1.15)
style(b, "(b) Proportions", "what sequencing reports")

# (c) apparent fold change of one unchanged taxon
c = axes[2]
i = 1                                           # any of the 29 unchanged taxa
prop_fc = (case[i] / case.sum()) / (ctrl[i] / ctrl.sum())
gm = lambda v: np.exp(np.log(v).mean())
clr_fc = (case[i] / gm(case)) / (ctrl[i] / gm(ctrl))
# the same calculation with 300 taxa (one bloomer, 299 unchanged)
o300 = np.exp(rng.normal(0, 0.8, 299)); o300 = 95 * o300 / o300.sum()
c300, k300 = np.r_[5.0, o300], np.r_[95.0, o300]
clr300 = (k300[1] / gm(k300)) / (c300[1] / gm(c300))
rows = [("proportions", prop_fc, CRIMSON), ("CLR, 30 taxa", clr_fc, DARK), ("CLR, 300 taxa", clr300, DARK),
        ("cells (total load\nmeasured)", 1.0, GREYF)]
Y = np.arange(len(rows))[::-1]
for y, (nm, v, col) in zip(Y, rows):
    c.barh(y, v, height=0.55, color=col)
    c.text(1.05, y, f"×{v:.2f}", ha="left", va="center", fontsize=9, fontweight="semibold",
           color=CRIMSON if v < 0.8 else TEXT)
c.axvline(1, color=TEXT, lw=0.8, ls=(0, (3, 2)))
c.set_yticks(Y); c.set_yticklabels([r[0] for r in rows], fontsize=8.8, color=TEXT)
c.set_xlim(0, 1.35); c.set_xticks([0, 0.5, 1]); c.tick_params(axis="y", length=0)
c.set_xlabel("apparent change (truth = ×1)", color=TEXT2, fontsize=9)
style(c, "(c) One unchanged taxon", "how big its 'change' looks")
fig.text(0.125, -0.06, "CLR compares each taxon with the average taxon (geometric mean), so one bloom moves the reference "
         "only a little — the more taxa, the less.", fontsize=8.8, color=TEXT2)
header(a, "Metagenomics · compositional data", "One taxon grew, and every other taxon seems to shrink", dy=40, x1=3.6)
fig.savefig(OUT + "compositional.png")
print("ok n5 prop", round(prop_fc, 3), "clr30", round(clr_fc, 3), "clr300", round(clr300, 3))
