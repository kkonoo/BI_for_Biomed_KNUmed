from gallery import *
from m_sim import missing_sim
from scipy import stats
import warnings; warnings.filterwarnings("ignore")

# Missing values: (a) which values go missing (known only in a simulation); (b) what filling them with a
# small constant does to a volcano plot
S = missing_sim(0)
X, O, tl, n = S["X"], S["O"], S["tl"], S["n"]
fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.5), gridspec_kw={"wspace": 0.28})
DARK = "#5E6670"; GREYF = "#B8BEC7"

a = axes[0]
bins = np.arange(15, 31, 0.4)
det = ~np.isnan(O)
h_det, _ = np.histogram(X[det], bins); h_low, _ = np.histogram(X[S["below"]], bins)
h_ch, _ = np.histogram(X[S["chance"]], bins)
ctr = (bins[:-1] + bins[1:]) / 2; w = 0.36
a.bar(ctr, h_det, width=w, color=GREYF, label="measured")
a.bar(ctr, h_low, width=w, bottom=h_det, color=CRIMSON, label="missing: below the detection limit (MNAR)")
a.bar(ctr, h_ch, width=w, bottom=h_det + h_low, color=DARK, label="missing: by chance (MAR)")
a.plot([S["lod"]] * 2, [0, (h_det + h_low + h_ch).max() * 1.15], color=CRIMSON, lw=0.9, ls=(0, (3, 2)))
a.text(S["lod"] - 0.2, (h_det + h_low + h_ch).max() * 1.08, "detection limit", ha="right", va="center", fontsize=8.5,
       color=CRIMSON)
a.set_ylim(0, (h_det + h_low + h_ch).max() * 1.62)
a.set_xlim(15.5, 30); a.set_yticks([])
for sp in ("left", "top", "right"): a.spines[sp].set_visible(False)
a.spines["bottom"].set_color(LINE); a.tick_params(colors=TEXT2, labelsize=8.3)
a.set_xlabel("true log₂ intensity", color=TEXT2, fontsize=9)
a.legend(frameon=False, fontsize=8.3, loc="upper left", labelcolor=TEXT, handlelength=1, bbox_to_anchor=(-0.02, 1.02))
miss = 100 * np.isnan(O).mean(); mnar = 100 * S["below"].sum() / np.isnan(O).sum()
a.set_title("(a) Which values go missing?", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
a.text(0, 1.03, f"simulated: {miss:.0f}% missing, {mnar:.0f}% of them because they were low", transform=a.transAxes,
       fontsize=9, color=TEXT2)

# (b) impute with a small constant, then t-test
b = axes[1]
const = np.nanmin(O)
I = np.where(np.isnan(O), const, O)
t, p = stats.ttest_ind(I[:, n:], I[:, :n], axis=1); p = np.where(np.isnan(p), 1, p)
lfc = I[:, n:].mean(1) - I[:, :n].mean(1)
imp = np.isnan(O).any(1)
lp = -np.log10(p)
b.axvspan(-0.8, 0.8, color=FAINT, lw=0, zorder=0)
YMAX = 7
lpc = np.minimum(lp, YMAX - 0.15)
b.scatter(lfc[~imp], lpc[~imp], s=6, color=GREYF, lw=0, alpha=0.8, zorder=2)
b.scatter(lfc[imp], lpc[imp], s=7, color=CRIMSON, lw=0, alpha=0.55, zorder=3)
for x_, y_ in zip(lfc[lp > YMAX], lp[lp > YMAX]):
    b.text(x_ - 0.1, YMAX - 0.15, "↑ p = 10⁻" + f"{y_:.0f}".translate(str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")), ha="right", va="center", fontsize=7.8, color=CRIMSON)
for v in (-1, 1): b.axvline(v, color=TEXT2, lw=0.7, ls=(0, (3, 2)))
b.axhline(-np.log10(0.05), color=TEXT2, lw=0.7, ls=(0, (3, 2)))
hit = (np.abs(lfc) > 1) & (p < 0.05)
b.text(0, 0.66, "every true\nchange is\ninside this\nband", transform=b.get_xaxis_transform(), ha="center", va="top",
       fontsize=8.3, color=TEXT2, linespacing=1.15)
b.text(0.98, 0.9, f"{hit.sum()} 'hits' (> 2-fold, p < 0.05)\n{(hit & imp).sum()} rest on imputed values\n"
       f"{(hit & (tl == 0)).sum()} have no true change at all",
       transform=b.transAxes, ha="right", va="top", fontsize=8.8, color=CRIMSON, fontweight="semibold",
       linespacing=1.25, bbox=dict(fc="white", ec="none", pad=1.5), zorder=5)
b.scatter([], [], s=14, color=CRIMSON, label="≥ 1 value imputed")
b.scatter([], [], s=14, color=GREYF, label="all measured")
b.legend(frameon=False, fontsize=8.3, loc="upper left", bbox_to_anchor=(-0.02, 1.0), labelcolor=TEXT,
         handletextpad=0.2)
b.set_xlabel("log₂ fold change (after imputation)", color=TEXT2, fontsize=9)
b.set_ylabel("−log₁₀ p", color=TEXT2, fontsize=9)
for sp in ("left", "bottom"): b.spines[sp].set_color(LINE)
b.tick_params(colors=TEXT2, labelsize=8.3)
b.set_ylim(0, YMAX)
b.set_title("(b) Fill the blanks with a small constant", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT,
            pad=22)
b.text(0, 1.03, f"same data; no protein truly changes more than 1.7-fold", transform=b.transAxes, fontsize=9, color=TEXT2)
header(a, "Proteomics · missing values", "A blank usually means 'too low to see', not 'random'", dy=40, x1=2.28)
fig.savefig(OUT + "missing_values.png")
print("ok m6", round(miss), round(mnar), "hits", hit.sum(), "hits with no true change", (hit & (tl == 0)).sum(),
      "imputed among hits", (hit & imp).sum(), "const", round(const, 1))
