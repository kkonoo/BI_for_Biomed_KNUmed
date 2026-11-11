from gallery import *
from n_sim import two_groups, bray_curtis, pcoa, permanova
rng = np.random.default_rng(111)

# PERMANOVA: p shrinks as samples are added, R² does not. One simulated pool (400 controls + 400 cases, a small
# shift in 40 of 150 taxa); subsamples of increasing size, 25 random subsamples per size (medians shown).
C, g = two_groups(400, rng, 0.8, n_eff=40)
D = bray_curtis(C)
F, R2_all, p_all = permanova(D, g, rng)
xy, expl = pcoa(D)

sizes = [10, 20, 40, 80, 160, 400]
med_p, med_r2 = [], []
for n in sizes:
    ps, rs = [], []
    for rep in range(25 if n < 400 else 1):
        idx = np.r_[rng.choice(400, n, replace=False), 400 + rng.choice(400, n, replace=False)]
        f, r2, p = permanova(D[np.ix_(idx, idx)], g[idx], rng)
        ps.append(p); rs.append(r2)
    med_p.append(np.median(ps)); med_r2.append(np.median(rs))

fig = plt.figure(figsize=(9.8, 4.6))
gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.15], height_ratios=[1, 1], wspace=0.3, hspace=0.35,
                      left=0.06, right=0.98, top=0.78, bottom=0.12)
DARK = "#5E6670"; GREYF = "#B8BEC7"
a = fig.add_subplot(gs[:, 0])
o = rng.permutation(len(g))
a.scatter(xy[o, 0], xy[o, 1], s=6, c=np.where(g[o] == 1, CRIMSON, GREYF), lw=0, alpha=0.75)
a.set_xticks([]); a.set_yticks([])
for sp in ("top", "right"): a.spines[sp].set_visible(False)
for sp in ("left", "bottom"): a.spines[sp].set_color(LINE)
a.set_xlabel(f"PCoA 1 ({expl[0]:.0%})", color=TEXT2, fontsize=9); a.set_ylabel(f"PCoA 2 ({expl[1]:.0%})", color=TEXT2,
                                                                                fontsize=9)
a.text(0.03, 0.97, f"PERMANOVA\np = 0.001\nR² = {R2_all:.3f}", transform=a.transAxes, ha="left", va="top",
       fontsize=9.5, color=CRIMSON, fontweight="semibold", linespacing=1.2)
a.set_title("(a) 800 people", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
a.text(0, 1.03, "grey = controls, red = cases (simulated)", transform=a.transAxes, fontsize=9, color=TEXT2)

N = 2 * np.array(sizes)
b1 = fig.add_subplot(gs[0, 1]); b2 = fig.add_subplot(gs[1, 1], sharex=b1)
b1.plot(N, med_p, color=DARK, lw=1.6, marker="o", ms=4)
b1.axhline(0.05, color=TEXT2, lw=0.7, ls=(0, (3, 2)))
b1.text(N[-1], 0.06, "0.05", ha="right", va="bottom", fontsize=8, color=TEXT2)
b1.set_yscale("log"); b1.set_ylim(5e-4, 1.2)
b1.set_yticks([0.001, 0.01, 0.1, 1]); b1.set_yticklabels(["0.001", "0.01", "0.1", "1"])
b1.set_ylabel("p", color=TEXT2, fontsize=9)
b2.plot(N, np.array(med_r2) * 100, color=CRIMSON, lw=1.6, marker="o", ms=4)
b2.plot(N, 100 / (N - 1), color=GREYF, lw=1.2, ls=(0, (3, 2)))
b2.text(N[-1], max(med_r2) * 100 * 1.2, "red = observed R²\ndashed = expected from noise alone", ha="right",
        va="top", fontsize=8.3, color=TEXT2, linespacing=1.2)
b2.axhline(R2_all * 100, color=CRIMSON, lw=0.6, alpha=0.4)
b2.set_ylabel("R² (%)", color=TEXT2, fontsize=9); b2.set_ylim(0, max(med_r2) * 100 * 1.25)
b2.set_xscale("log"); b2.set_xticks(N); b2.set_xticklabels([str(n) for n in N])
b2.set_xlabel("people in the study", color=TEXT2, fontsize=9)
for ax in (b1, b2):
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"): ax.spines[sp].set_color(LINE)
    ax.tick_params(colors=TEXT2, labelsize=8.3); ax.minorticks_off()
plt.setp(b1.get_xticklabels(), visible=False)
b1.set_title("(b) Add people", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
b1.text(0, 1.06, "median of 25 random subsets per size", transform=b1.transAxes, fontsize=9, color=TEXT2)
header(a, "Metagenomics · PERMANOVA", "With enough people p is tiny; R² says how much the grouping explains",
       dy=40, x1=2.3)
fig.savefig(OUT + "permanova.png")
print("ok n11 R2_all", round(R2_all, 4), "p", p_all, "med_p", np.round(med_p, 3), "med_r2", np.round(med_r2, 3))
