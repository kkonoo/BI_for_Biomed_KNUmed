from gallery import *
rng = np.random.default_rng(4)

# On the log scale the two group means are β0 and β0 + β1; β1 = log2 fold change; Wald test z = β1 / SE
fig = plt.figure(figsize=(9.0, 3.9))
ax = fig.add_axes([0.07, 0.13, 0.48, 0.6])
for sp in ("left", "bottom"): ax.spines[sp].set_color(LINE)
ax.tick_params(colors=TEXT2, labelsize=8.5)
GREYF = "#9AA1AB"
b0, b1 = 6.0, 1.6
for x0, m, col in ((0, b0, GREYF), (1, b0 + b1, CRIMSON)):
    d = rng.normal(0, 0.28, 4); y = m + d - d.mean()
    ax.scatter(x0 + np.array([-0.06, -0.02, 0.02, 0.06]), y, s=48, color=col, edgecolor="white", lw=0.8, zorder=3)
    ax.plot([x0 - 0.2, x0 + 0.2], [m, m], color=TEXT, lw=1.6)
ax.plot([0.2, 1.45], [b0, b0], color=LINE, lw=0.9, ls=(0, (3, 2)))
ax.annotate("", (1.42, b0 + b1), (1.42, b0), arrowprops=dict(arrowstyle="<->", color=CRIMSON, lw=1.3))
ax.errorbar([1.42], [b0 + b1], yerr=[[0.35], [0.35]], fmt="none", ecolor=tint(0.5), elinewidth=6, capsize=0,
            zorder=1)
ax.text(1.5, b0 + b1 / 2, "β₁ = log₂FC\n= 1.6 (≈ 3×)", ha="left", va="center", fontsize=9.3, color=CRIMSON,
        fontweight="semibold", linespacing=1.25)
ax.text(1.5, b0 + b1 + 0.35, "± SE", ha="left", va="center", fontsize=8.4, color=TEXT2)
ax.text(-0.27, b0, "β₀", ha="right", va="center", fontsize=11, fontweight="bold", color=TEXT)
ax.set_xlim(-0.45, 2.05); ax.set_ylim(4.9, 8.4)
ax.set_xticks([0, 1]); ax.set_xticklabels(["control", "case"]); ax.set_yticks([])
ax.set_ylabel("log₂ expression", color=TEXT2, fontsize=9)

tx = fig.add_axes([0.6, 0.13, 0.38, 0.6]); tx.axis("off"); tx.set_xlim(0, 1); tx.set_ylim(0, 1)
tx.text(0, 0.9, "log₂(μ) = β₀ + β₁·condition", fontsize=12.5, fontweight="bold", color=TEXT, va="center")
tx.text(0, 0.74, "β₁ = log₂(μ_case / μ_control)", fontsize=9.6, color=TEXT, va="center")
tx.plot([0, 0.95], [0.62, 0.62], color=LINE, lw=0.9)
tx.text(0, 0.5, "Wald test", fontsize=10.5, fontweight="bold", color=CRIMSON, va="center")
tx.text(0, 0.37, "z = β₁ / SE(β₁) → p-value", fontsize=9.6, color=TEXT, va="center")
tx.text(0, 0.21, "LRT", fontsize=10.5, fontweight="bold", color=TEXT, va="center")
tx.text(0, 0.08, "full vs reduced model → 3+ groups, time course", fontsize=9.2, color=TEXT, va="center")
header(ax, "Differential expression · the test", "The log₂ fold change is a coefficient; the test asks if it is 0",
       dy=8, x1=2.05)
fig.savefig(OUT + "glm_lfc.png"); print("ok h26 glm lfc")
