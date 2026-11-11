from gallery import *
from m_sim import missing_sim
import warnings; warnings.filterwarnings("ignore")

# Three ways to handle missing values on the same simulated data (as in missing_values.png):
# estimated vs true fold change for each protein
S = missing_sim(0)
O, tl, n, rng = S["O"], S["tl"], S["n"], S["rng"]
obs = ~np.isnan(O)
I_min = np.where(np.isnan(O), np.nanmin(O), O)
I_ds = O.copy()                                             # Perseus-style: draw from a down-shifted normal per sample
for j in range(2 * n):
    col = O[:, j]; s = np.nanstd(col); m = np.nanmean(col); na = np.isnan(col)
    I_ds[na, j] = rng.normal(m - 1.8 * s, 0.3 * s, na.sum())
keep_f = (obs[:, :n].sum(1) >= 4) & (obs[:, n:].sum(1) >= 4)
allk = np.ones(len(tl), bool)
cases = [("(a) Fill with the minimum", "a small constant", I_min, allk),
         ("(b) Draw low random values", "down-shifted normal (Perseus)", I_ds, allk),
         ("(c) Filter, then no imputation", "≥ 4 of 6 measured in each group", O, keep_f)]
fig, axes = plt.subplots(1, 3, figsize=(9.8, 4.2), gridspec_kw={"wspace": 0.12})
DARK = "#5E6670"; GREYF = "#B8BEC7"
jit = rng.uniform(-0.04, 0.04, len(tl))
for ax, (t, sub, I, keep) in zip(axes, cases):
    est = np.nanmean(I[keep][:, n:], 1) - np.nanmean(I[keep][:, :n], 1)
    tru = tl[keep]; bad = np.abs(est - tru) > 1
    ax.plot([-1, 1], [-1, 1], color=TEXT2, lw=0.8, ls=(0, (3, 2)))
    ax.scatter(tru[~bad] + jit[keep][~bad], est[~bad], s=5, color=GREYF, lw=0, alpha=0.6)
    ax.scatter(tru[bad] + jit[keep][bad], est[bad], s=6, color=CRIMSON, lw=0, alpha=0.6)
    ax.set_xlim(-1.1, 1.1); ax.set_ylim(-4.6, 4.6)
    for sp in ("left", "bottom"): ax.spines[sp].set_color(LINE)
    ax.tick_params(colors=TEXT2, labelsize=8.3)
    ax.set_xlabel("true log₂ fold change", color=TEXT2, fontsize=9)
    if ax is axes[0]: ax.set_ylabel("estimated log₂ fold change", color=TEXT2, fontsize=9)
    else: ax.set_yticklabels([])
    ax.text(0.04, 0.97, f"{keep.sum():,} proteins kept", transform=ax.transAxes, ha="left", va="top", fontsize=9,
            color=TEXT)
    ax.text(0.04, 0.89, f"{bad.sum()} off by > 2-fold", transform=ax.transAxes, ha="left", va="top", fontsize=9.5,
            fontweight="bold", color=CRIMSON if bad.sum() else TEXT)
    ax.set_title(t, loc="left", fontsize=11.8, fontweight="semibold", color=TEXT, pad=22)
    ax.text(0, 1.03, sub, transform=ax.transAxes, fontsize=8.8, color=TEXT2)
fig.text(0.125, -0.04, "Filtering keeps fewer proteins, but their fold changes can be trusted. Whatever you choose: "
         "report it, and check that the top hits survive without imputation.", fontsize=8.8, color=TEXT2)
header(axes[0], "Proteomics · imputation", "How you fill the blanks decides the fold changes you report", dy=40,
       x1=3.25)
fig.savefig(OUT + "imputation.png")
print("ok m12", [(int(k.sum())) for *_, k in cases])
