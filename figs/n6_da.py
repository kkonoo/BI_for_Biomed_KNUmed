from gallery import *
from n_da import simulate, methods
import warnings; warnings.filterwarnings("ignore")
rng = np.random.default_rng(2)

# Differential abundance: five simplified approaches on one simulated dataset (40 + 40 people, 120 taxa,
# 13 truly changed in cell numbers), plus a reference that knows the total load.
S = simulate(rng); M = methods(S, rng); truth = S["truth"]
names = list(M.keys())
calls = np.array([M[k] for k in names])
test = calls[:-1]                                   # the five approaches (reference row excluded from agreement)
show = np.where(truth | calls.any(0))[0]
order = sorted(show, key=lambda j: (not truth[j], -calls[:, j].sum()))
fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.4), gridspec_kw={"wspace": 0.08, "width_ratios": [1.7, 1]})
DARK = "#5E6670"; GREYF = "#B8BEC7"

a = axes[0]
nr = len(names)
for r in range(nr):
    y = nr - 1 - r
    for cidx, j in enumerate(order):
        c = calls[r, j]
        col = (DARK if truth[j] else CRIMSON) if c else FAINT
        a.add_patch(mp.Rectangle((cidx, y + 0.1), 0.86, 0.8, fc=col, ec="none"))
ntrue = truth.sum()
a.plot([0, ntrue - 0.14], [nr + 0.25] * 2, color=DARK, lw=2)
a.text(0, nr + 0.42, f"{ntrue} taxa truly changed", ha="left", va="bottom", fontsize=8.5, color=DARK,
       fontweight="semibold")
a.plot([ntrue, len(order) - 0.14], [nr + 0.25] * 2, color=CRIMSON, lw=2)
a.text(ntrue, nr + 0.42, "unchanged, yet called by at least one approach", ha="left", va="bottom", fontsize=8.5,
       color=CRIMSON, fontweight="semibold")
a.plot([-0.2, len(order)], [1.0, 1.0], color=TEXT2, lw=0.6, ls=(0, (3, 2)))
a.set_xlim(-0.3, len(order) + 0.2); a.set_ylim(-0.1, nr + 1.3)
a.set_yticks(np.arange(nr) + 0.5); a.set_yticklabels(names[::-1], fontsize=8.4, color=TEXT, linespacing=1.05)
a.set_xticks([]); a.tick_params(axis="y", length=0)
for sp in a.spines.values(): sp.set_visible(False)
a.set_xlabel("one column = one taxon; filled = called significant", color=TEXT2, fontsize=8.8)
a.set_title("(a) Who calls what", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
a.text(0, 1.03, "simplified versions of each approach, not the packages", transform=a.transAxes, fontsize=9,
       color=TEXT2)

b = axes[1]
tp = (calls & truth).sum(1); fp = (calls & ~truth).sum(1)
Y = np.arange(nr)[::-1] + 0.5
b.barh(Y, tp, height=0.6, color=DARK)
b.barh(Y, fp, left=tp, height=0.6, color=CRIMSON)
for y, t_, f_ in zip(Y, tp, fp):
    b.text(t_ + f_ + 0.6, y, f"{t_} + {f_}", ha="left", va="center", fontsize=8.6, color=TEXT)
b.set_ylim(-0.1, nr + 1.3); b.set_yticks([])
b.set_xlim(0, (tp + fp).max() * 1.25)
for sp in ("top", "right", "left"): b.spines[sp].set_visible(False)
b.spines["bottom"].set_color(LINE); b.tick_params(colors=TEXT2, labelsize=8.3)
b.set_xlabel("taxa called: right + wrong", color=TEXT2, fontsize=8.8)
allc = (test.sum(0) == len(test)).sum(); one = (test.sum(0) == 1).sum()
b.text(0.98, 0.98, f"called by all five: {allc}\ncalled by only one: {one}", transform=b.transAxes, ha="right",
       va="top", fontsize=8.8, color=TEXT, linespacing=1.3)
b.set_title("(b) Right and wrong calls", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
b.text(0, 1.03, "dark = truly changed, red = unchanged", transform=b.transAxes, fontsize=9, color=TEXT2)
header(a, "Metagenomics · differential abundance", "Same data, different method, different answer", dy=40, x1=1.62)
fig.savefig(OUT + "da_disagreement.png")
print("ok n6", dict(zip([n.replace("\n", " ") for n in names], zip(tp, fp))), "all", allc, "one", one)
