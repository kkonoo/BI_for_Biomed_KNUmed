from gallery import *
rng = np.random.default_rng(8)

# Counts are overdispersed (left); dispersion is shrunk toward a trend fitted from all genes (right)
fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.6), gridspec_kw={"wspace": 0.22})
DARK = "#5E6670"; GREYF = "#B8BEC7"
G = 4000
mu = np.exp(rng.normal(3.4, 1.9, G))
disp = 0.02 + 18 / (mu + 30)
var = mu + disp * mu ** 2
obs = var * np.exp(rng.normal(0, 0.28, G))

def style(ax, title, sub):
    for sp in ("left", "bottom"): ax.spines[sp].set_color(LINE)
    ax.tick_params(colors=TEXT2, labelsize=8.5)
    ax.set_title(title, loc="left", fontsize=12.5, fontweight="semibold", color=TEXT, pad=24)
    ax.text(0, 1.03, sub, transform=ax.transAxes, ha="left", va="bottom", fontsize=9.3, color=TEXT2)

# (a) overdispersion
ax = axes[0]
ax.scatter(mu, obs, s=4, color=GREYF, alpha=0.45, lw=0)
m = np.logspace(0, 5, 200)
ax.plot(m, m, color=DARK, lw=1.6, ls=(0, (4, 2)))
ax.plot(m, m + (0.02 + 18 / (m + 30)) * m ** 2, color=CRIMSON, lw=2)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(1, 1e5); ax.set_ylim(0.3, 3e8)
ax.set_xlabel("mean count across replicates", color=TEXT2, fontsize=9)
ax.set_ylabel("variance across replicates", color=TEXT2, fontsize=9)
ax.text(5, 4e6, "negative binomial\nvar = μ + φμ²", ha="left", va="center", fontsize=9.3,
        color=CRIMSON, fontweight="semibold", linespacing=1.3)
ax.text(3e4, 6e3, "Poisson\nvar = μ", ha="left", va="center", fontsize=9.3, color=DARK, fontweight="semibold",
        linespacing=1.3)
ax.text(9e4, 0.6, "replicates vary far more than Poisson\n→ t-test or Poisson on counts\n   = false positives",
        ha="right", va="bottom", fontsize=8.8, color=TEXT2, linespacing=1.4)
style(ax, "(a) Counts are overdispersed", "each dot = one gene")

# (b) shrinkage
ax = axes[1]
raw = disp * np.exp(rng.normal(0, 0.55, G))
trend = 0.02 + 18 / (mu + 30)
shrunk = np.exp(0.35 * np.log(raw) + 0.65 * np.log(trend))
o = np.argsort(mu)
ax.scatter(mu, raw, s=4, color=GREYF, alpha=0.45, lw=0)
ax.scatter(mu, shrunk, s=4, color=tint(0.55), alpha=0.6, lw=0)
ax.plot(mu[o], trend[o], color=CRIMSON, lw=2.2)
# one gene's estimate pulled toward the trend
k = np.argmin(np.abs(np.log(mu) - np.log(40)) + 5 * (raw < trend * 2.6))
ax.annotate("", (mu[k], shrunk[k]), (mu[k], raw[k]), arrowprops=dict(arrowstyle="-|>", color=TEXT, lw=1.3,
                                                                    mutation_scale=10))
ax.plot([mu[k]], [raw[k]], "o", ms=6, mfc="white", mec=TEXT, mew=1.3, zorder=5)
ax.plot([mu[k]], [shrunk[k]], "o", ms=6, color=CRIMSON, mec="white", mew=1, zorder=5)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(1, 1e5); ax.set_ylim(3e-3, 30)
ax.set_xlabel("mean count", color=TEXT2, fontsize=9)
ax.set_ylabel("dispersion φ", color=TEXT2, fontsize=9)
ax.text(1.4, 0.0065, "grey = one gene alone (n = 3: noisy)\npink = pulled toward the trend (used)\nline = trend fitted from all genes",
        ha="left", va="bottom", fontsize=8.8, color=TEXT2, linespacing=1.45)
ax.text(mu[k] * 1.35, raw[k], "one gene", ha="left", va="center", fontsize=8.8, color=TEXT)
style(ax, "(b) Borrowing from all genes", "why DESeq2 / edgeR work with n = 3")

header(axes[0], "Differential expression · count model",
       "Counts vary more than Poisson, and n = 3 is rescued by shrinkage", dy=44, x1=2.22)
fig.savefig(OUT + "mean_variance.png"); print("ok h3")
