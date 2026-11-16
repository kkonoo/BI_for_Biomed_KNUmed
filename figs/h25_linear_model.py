from gallery import *
from scipy.stats import norm

# Two figures: (1) a linear model = mean part + noise part; (2) variance decomposition with and without batch
GREYF = "#9AA1AB"; DARK = "#5E6670"; rng = np.random.default_rng(2)

# 1 — linear model with a 0/1 condition
fig = plt.figure(figsize=(9.0, 3.9))
ax = fig.add_axes([0.07, 0.13, 0.5, 0.6])
for sp in ("left", "bottom"): ax.spines[sp].set_color(LINE)
ax.tick_params(colors=TEXT2, labelsize=8.5)
b0, b1, s = 4.0, 2.2, 0.55
for x0, col in ((0, GREYF), (1, CRIMSON)):
    y = b0 + b1 * x0 + rng.normal(0, s, 7)
    ax.scatter(x0 + rng.uniform(-0.04, 0.04, 7), y, s=40, color=col, edgecolor="white", lw=0.7, zorder=3)
    yy = np.linspace(b0 + b1 * x0 - 3 * s, b0 + b1 * x0 + 3 * s, 100)
    ax.fill_betweenx(yy, x0 + 0.12, x0 + 0.12 + norm.pdf(yy, b0 + b1 * x0, s) * 0.32, color=tint(0.18), lw=0)
    ax.plot(x0 + 0.12 + norm.pdf(yy, b0 + b1 * x0, s) * 0.32, yy, color=tint(0.6), lw=1)
ax.plot([-0.2, 1.2], [b0 - 0.2 * b1, b0 + 1.2 * b1], color=TEXT, lw=1.6, zorder=2)
ax.annotate("", (1.42, b0 + b1), (1.42, b0), arrowprops=dict(arrowstyle="<->", color=CRIMSON, lw=1.2))
ax.plot([0, 1.42], [b0, b0], color=LINE, lw=0.9, ls=(0, (3, 2)))
ax.text(1.47, b0 + b1 / 2, "β₁", ha="left", va="center", fontsize=11, color=CRIMSON, fontweight="bold")
ax.text(-0.27, b0, "β₀", ha="right", va="center", fontsize=11, color=TEXT, fontweight="bold")
ax.text(0.33, b0 - 1.75, "same spread σ²\nin both groups", ha="left", va="center", fontsize=8.6, color=TEXT2,
        linespacing=1.25)
ax.set_xlim(-0.45, 1.75); ax.set_ylim(1.8, 8.4)
ax.set_xticks([0, 1]); ax.set_xticklabels(["control (x = 0)", "case (x = 1)"])
ax.set_ylabel("expression", color=TEXT2, fontsize=9); ax.set_yticks([])
tx = fig.add_axes([0.62, 0.13, 0.36, 0.64]); tx.axis("off"); tx.set_xlim(0, 1); tx.set_ylim(0, 1)
tx.text(0, 0.88, "y  =  β₀ + β₁·x  +  ε", fontsize=14, fontweight="bold", color=TEXT, va="center")
tx.text(0, 0.72, "mean part = β₀ + β₁·x   ·   noise part = ε", fontsize=9, color=TEXT2, va="center")
rows = (("β₀", "mean of the control group"), ("β₁", "condition effect = what we test"),
        ("ε", "noise, spread σ² (the same everywhere)"))
for i, (k, v) in enumerate(rows):
    tx.text(0, 0.5 - i * 0.16, k, fontsize=11, fontweight="bold", color=CRIMSON if k == "β₁" else TEXT, va="center")
    tx.text(0.1, 0.5 - i * 0.16, v, fontsize=9.2, color=TEXT, va="center")
header(ax, "Differential expression · regression", "A linear model has a mean part and a noise part", dy=8,
       x1=1.88)
fig.savefig(OUT + "linear_model.png"); plt.close(fig)

# 2 — variance decomposition
fig = plt.figure(figsize=(9.0, 3.6))
ax = fig.add_axes([0.2, 0.12, 0.62, 0.62]); ax.axis("off"); ax.set_xlim(0, 100); ax.set_ylim(-0.4, 2.2)
bars = (("~ condition", [("condition", 25, CRIMSON), ("residual (batch hidden inside)", 75, "#D5D9DF")]),
        ("~ batch + condition", [("condition", 25, CRIMSON), ("batch", 45, tint(0.35)), ("residual", 30, "#D5D9DF")]))
for i, (lab, parts) in enumerate(bars):
    y = 1.5 - i * 1.15; left = 0
    ax.text(-2, y, lab, ha="right", va="center", fontsize=10, color=TEXT, family="DejaVu Sans Mono")
    for name, v, col in parts:
        ax.barh(y, v, left=left, height=0.5, color=col, lw=0)
        ax.text(left + v / 2, y, f"{name}\n{v}%", ha="center", va="center", fontsize=8.6, linespacing=1.2,
                color="white" if col in (CRIMSON,) else TEXT, fontweight="semibold")
        left += v
    ratio = parts[0][1] / parts[-1][1]
    ax.text(103, y, f"signal : noise\n= 1 : {1 / ratio:.1f}", ha="left", va="center", fontsize=9,
            color=CRIMSON if i else TEXT2, fontweight="semibold" if i else "normal", linespacing=1.25)
ax.text(0, -0.25, "Total variance = condition + batch + residual   ·   the same data, only the model differs",
        ha="left", va="center", fontsize=9, color=TEXT2)
header(ax, "Differential expression · variance decomposition", "Name the batch, and the condition stands out",
       x=-0.27, x1=1.25)
fig.savefig(OUT + "variance_decomp.png"); print("ok h25 linear model + variance decomposition")
