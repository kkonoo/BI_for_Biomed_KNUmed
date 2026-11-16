from gallery import *
rng = np.random.default_rng(5)

# PCA on VST data: what you hope to see, and what a batch effect looks like
fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.3), gridspec_kw={"wspace": 0.22})
GREYF = "#9AA1AB"

def style(ax, title, sub, xl, yl):
    for sp in ("left", "bottom"): ax.spines[sp].set_color(LINE)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel(xl, color=TEXT2, fontsize=9); ax.set_ylabel(yl, color=TEXT2, fontsize=9)
    ax.set_title(title, loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=24)
    ax.text(0, 1.03, sub, transform=ax.transAxes, ha="left", va="bottom", fontsize=9.1, color=TEXT2)
    ax.set_xlim(-6, 6); ax.set_ylim(-4.2, 4.2)

def pts(ax, cx, cy, cond, batch, n=3, sd=0.45):
    x = cx + rng.normal(0, sd, n); y = cy + rng.normal(0, sd, n)
    col = CRIMSON if cond == "case" else GREYF
    if batch == 1:
        ax.scatter(x, y, s=70, color=col, edgecolor="white", lw=0.8, zorder=3)
    else:
        ax.scatter(x, y, s=70, facecolor="white", edgecolor=col, lw=1.8, zorder=3)

# (a) condition drives PC1, batch is minor
ax = axes[0]
pts(ax, -3.2, 0.5, "control", 1); pts(ax, -3.0, -0.6, "control", 2)
pts(ax, 3.1, 0.6, "case", 1); pts(ax, 3.3, -0.5, "case", 2)
ax.text(-3.1, 2.4, "control", ha="center", fontsize=9.5, color=TEXT2, fontweight="semibold")
ax.text(3.2, 2.4, "case", ha="center", fontsize=9.5, color=CRIMSON, fontweight="semibold")
style(ax, "(a) What you hope to see", "groups split on PC1, replicates close together",
      "PC1 (58% of variance)", "PC2 (9%)")

# (b) batch drives PC1, condition only on PC2
ax = axes[1]
pts(ax, -3.2, -1.5, "control", 1); pts(ax, -3.0, 1.6, "case", 1)
pts(ax, 3.0, -1.4, "control", 2); pts(ax, 3.3, 1.5, "case", 2)
ax.plot([0, 0], [-3.8, 3.8], color=LINE, lw=1, ls=(0, (3, 3)))
ax.text(-3.1, 3.3, "batch 1", ha="center", fontsize=9.5, color=TEXT, fontweight="semibold")
ax.text(3.15, 3.3, "batch 2", ha="center", fontsize=9.5, color=TEXT, fontweight="semibold")
style(ax, "(b) Batch effect", "samples split by processing day, not by biology",
      "PC1 (41%) = batch", "PC2 (14%) = condition")

# legend
lx, ly = 0.0, -0.2
axes[0].scatter([lx], [ly], transform=axes[0].transAxes, s=55, color=CRIMSON, edgecolor="white", clip_on=False)
axes[0].text(lx + 0.03, ly, "case", transform=axes[0].transAxes, va="center", fontsize=8.6, color=TEXT2)
axes[0].scatter([lx + 0.16], [ly], transform=axes[0].transAxes, s=55, color=GREYF, edgecolor="white", clip_on=False)
axes[0].text(lx + 0.19, ly, "control", transform=axes[0].transAxes, va="center", fontsize=8.6, color=TEXT2)
axes[0].scatter([lx + 0.36], [ly], transform=axes[0].transAxes, s=55, color=TEXT2, edgecolor="white", clip_on=False)
axes[0].text(lx + 0.39, ly, "filled = batch 1", transform=axes[0].transAxes, va="center", fontsize=8.6, color=TEXT2)
axes[0].scatter([lx + 0.66], [ly], transform=axes[0].transAxes, s=55, facecolor="white", edgecolor=TEXT2, lw=1.6,
                clip_on=False)
axes[0].text(lx + 0.69, ly, "open = batch 2", transform=axes[0].transAxes, va="center", fontsize=8.6, color=TEXT2)
axes[1].text(1.0, -0.2, "→ put batch in the design formula (§4.4)", transform=axes[1].transAxes, ha="right",
             va="center", fontsize=9, color=CRIMSON, fontweight="semibold")

header(axes[0], "Pre-processing · EDA", "Look at the samples before you test the genes", dy=44, x1=2.22)
fig.savefig(OUT + "eda_pca.png"); print("ok h23 pca")
