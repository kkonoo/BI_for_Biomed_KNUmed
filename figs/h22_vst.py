from gallery import *
rng = np.random.default_rng(11)

# log2(count + 1) leaves low-count genes noisy; DESeq2's VST makes the spread roughly equal at every level
G, N = 6000, 8
mu = np.exp(rng.normal(3.0, 2.0, G)).clip(0.3, 2e5)
a0, a1 = 0.05, 1.0                                    # dispersion trend α(μ) = a0 + a1/μ
alpha = a0 + a1 / mu
counts = rng.negative_binomial(1 / alpha[:, None], 1 / (1 + alpha[:, None] * mu[:, None]), size=(G, N))
logc = np.log2(counts + 1)
vst = np.log((1 + a1 + 2 * a0 * counts + 2 * np.sqrt(a0 * counts * (1 + a1 + a0 * counts))) / (4 * a0)) / np.log(2)
keep = counts.sum(1) > 0
GREYF = "#B8BEC7"

fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.0), sharey=True, gridspec_kw={"wspace": 0.12})
def style(ax, title, sub):
    for sp in ("left", "bottom"): ax.spines[sp].set_color(LINE)
    ax.tick_params(colors=TEXT2, labelsize=8.5)
    ax.set_title(title, loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=24)
    ax.text(0, 1.03, sub, transform=ax.transAxes, ha="left", va="bottom", fontsize=9.1, color=TEXT2)

for ax, X, title, sub, col in ((axes[0], logc, "(a) log₂(count + 1)", "low-count genes look the most variable", TEXT2),
                               (axes[1], vst, "(b) VST (DESeq2)", "spread roughly the same at every level", CRIMSON)):
    m = X[keep].mean(1); sd = X[keep].std(1, ddof=1)
    r = np.argsort(np.argsort(m)) / keep.sum() * 100
    ax.scatter(r, sd, s=3, color=GREYF, alpha=0.35, lw=0)
    o = np.argsort(r); bins = np.linspace(0, 100, 41); c = (bins[1:] + bins[:-1]) / 2
    med = [np.median(sd[(r >= lo) & (r < hi)]) for lo, hi in zip(bins[:-1], bins[1:])]
    ax.plot(c, med, color=CRIMSON if col == CRIMSON else "#5E6670", lw=2.4)
    ax.set_xlim(0, 100); ax.set_ylim(0, 1.85)
    ax.set_xlabel("genes ranked by mean count (low → high)", color=TEXT2, fontsize=9)
    style(ax, title, sub)
axes[0].set_ylabel("SD across samples", color=TEXT2, fontsize=9)
axes[0].text(3, 1.68, "noise at low counts\ndominates PCA and clustering", ha="left", va="center", fontsize=8.8,
             color=TEXT, linespacing=1.35)
axes[1].text(3, 1.68, "every gene contributes on a\ncomparable scale", ha="left", va="center", fontsize=8.8,
             color=CRIMSON, fontweight="semibold", linespacing=1.35)

header(axes[0], "Pre-processing · EDA", "Transform before you look: VST or rlog, for plots only", dy=44, x1=2.12)
fig.savefig(OUT + "vst_meansd.png"); print("ok h22 vst")
