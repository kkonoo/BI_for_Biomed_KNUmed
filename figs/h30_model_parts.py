from gallery import *

# A regression model = a mean model (where is each group's average? → β, β₁ = the effect) + a variance model
# (how far do values scatter around it? → σ², or μ + αμ² in DESeq2). Variance decomposition splits the total
# variation into the parts the mean model explains (condition, batch) and the residual the variance model takes.
# Two figures with the same layout: (1) linear regression, (2) the negative binomial GLM of DESeq2.
DARK = "#5E6670"; GREYF = "#EEF0F3"; GREYE = "#AEB4BD"; ARROW = "#9AA1AB"; MEAN_E = tint(0.45)
W, TOP = 10.0, 1.05


def box(ax, xc, yc, w, h, fc, ec, lw=1.0, r=0.1):
    ax.add_patch(mp.FancyBboxPatch((xc - w / 2, yc - h / 2), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                                   fc=fc, ec=ec, lw=lw, zorder=2))


def arrow(ax, x, y0, y1, c=ARROW):
    ax.annotate("", (x, y1), (x, y0), arrowprops=dict(arrowstyle="-|>", color=c, lw=1.3, mutation_scale=13,
                                                      shrinkA=0, shrinkB=0), zorder=1)


def dash(ax, xs, ys, c):
    ax.plot(xs, ys, color=c, lw=1.1, ls=(0, (3, 2.5)), zorder=1)


def figure(name, eyebrow, title, top, mean_eq, var_eq, mean_leaves, var_leaves, residual, note=None):
    h = 7.9 + (0.95 if note else 0); H = h + TOP
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0, 0, 1, h / H]); ax.set_xlim(0, W); ax.set_ylim(0, h); ax.axis("off")
    L, R = 2.45, 7.55                                   # centres of the two columns
    # top: the whole model
    y = h - 0.45
    box(ax, 5.0, y, 6.0, 0.66, FAINT, LINE)
    ax.text(5.0, y, top, ha="center", va="center", fontsize=16, fontweight="bold", color=TEXT, zorder=3)
    for x, q in ((L, "Where is the mean?"), (R, "How spread out?")):
        arrow(ax, x, y - 0.33, y - 1.08)
        ax.text(x + 0.14, y - 0.72, q, ha="left", va="center", fontsize=11.3, color=TEXT2)
    # ① mean model | ② variance model
    y -= 1.55
    for x, fc, t, s in ((L, CRIMSON, "①  Mean model", "the expected value in each group"),
                        (R, DARK, "②  Variance model", "how far values scatter around μ")):
        box(ax, x, y, 4.7, 0.9, fc, "none")
        ax.text(x, y + 0.17, t, ha="center", va="center", fontsize=14.5, fontweight="bold", color="white", zorder=3)
        ax.text(x, y - 0.2, s, ha="center", va="center", fontsize=11.3, color="white", alpha=0.92, zorder=3)
    y -= 0.95
    box(ax, L, y, 4.7, 0.66, tint(0.08), MEAN_E)
    box(ax, R, y, 4.7, 0.66, GREYF, GREYE)
    ax.text(L, y, mean_eq, ha="center", va="center", fontsize=13.8, fontweight="bold", color=CRIMSON, zorder=3)
    ax.text(R, y, var_eq, ha="center", va="center", fontsize=13.8, fontweight="bold", color=DARK, zorder=3)
    # leaves: coefficients | variance terms
    ytop, y = y - 0.33, y - 1.35
    yl = y - 0.62                                      # label line under the leaves
    for x, sym, lab, key in mean_leaves:
        arrow(ax, x, ytop, y + 0.31, c=MEAN_E)
        box(ax, x, y, 1.25, 0.62, CRIMSON if key else tint(0.08), "none" if key else MEAN_E)
        ax.text(x, y, sym, ha="center", va="center", fontsize=15.5, fontweight="bold",
                color="white" if key else CRIMSON, zorder=3)
        ax.text(x, yl, lab, ha="center", va="center", fontsize=11.3, linespacing=1.2,
                color=CRIMSON if key else TEXT2, fontweight="bold" if key else "normal")
    for x, sym, lab in var_leaves:
        arrow(ax, x, ytop, y + 0.31, c=GREYE)
        box(ax, x, y, 1.25, 0.62, GREYF, GREYE)
        ax.text(x, y, sym, ha="center", va="center", fontsize=15.5, fontweight="bold", color=DARK, zorder=3)
        ax.text(x, yl, lab, ha="center", va="center", fontsize=11.3, color=TEXT2, linespacing=1.2)
    # variance decomposition
    yd = yl - 1.05
    for x, sym, lab, key in mean_leaves[1:]:
        dash(ax, [x, x], [yl - 0.3, yd + 0.31], MEAN_E)
    vx = [x for x, *_ in var_leaves]
    if len(vx) == 1:
        dash(ax, [vx[0], vx[0]], [yl - 0.3, yd + 0.31], GREYE)
    else:
        yb = yl - 0.55
        for x in vx:
            dash(ax, [x, x], [yl - 0.3, yb], GREYE)
        dash(ax, [vx[0], vx[-1]], [yb, yb], GREYE)
        dash(ax, [R, R], [yb, yd + 0.31], GREYE)
    box(ax, 5.0, yd, 9.8, 0.62, FAINT, LINE)
    ax.text(5.0, yd, "Variance decomposition:   total  =  condition  +  batch  +  residual", ha="center",
            va="center", fontsize=13.8, fontweight="bold", color=TEXT, zorder=3)
    # where each part goes
    yb = yd - 1.25
    for x, t, s, fc, ec, c in ((1.7, "Condition", "variance explained by β₁", tint(0.08), MEAN_E, CRIMSON),
                               (5.0, "Batch", "variance explained by β₂", tint(0.08), MEAN_E, CRIMSON),
                               (8.3, "Residual", residual, GREYF, GREYE, DARK)):
        arrow(ax, x, yd - 0.31, yb + 0.48)
        box(ax, x, yb, 3.1, 0.96, fc, ec)
        ax.text(x, yb + 0.2, t, ha="center", va="center", fontsize=13.8, fontweight="bold", color=c, zorder=3)
        ax.text(x, yb - 0.2, s, ha="center", va="center", fontsize=11.3, color=TEXT2, zorder=3)
    if note:
        ax.text(5.0, yb - 1.02, note, ha="center", va="center", fontsize=11.3, color=CRIMSON, fontweight="semibold",
                linespacing=1.35)
    header(ax, eyebrow, title, dy=6)
    fig.savefig(OUT + name)


# 1 — linear regression
figure("regression_parts.png", "Differential expression · regression",
       "A regression model has two parts: the mean and the spread",
       "Y  =  μ + ε,      ε ~ Normal(0, σ²)",
       "μ = β₀ + β₁·condition + β₂·batch", "Var = σ²  (one value)",
       [(0.95, "β₀", "control mean", False), (2.45, "β₁", "= condition effect", True),
        (3.95, "β₂", "batch effect", False)],
       [(7.55, "σ²", "the same for every sample")],
       "what is left: ε ~ N(0, σ²)")

# 2 — DESeq2: negative binomial GLM
figure("nbglm_parts.png", "Differential expression · DESeq2",
       "DESeq2 keeps the same two parts, fitted to counts",
       "Y  ~  NegBinomial(μ, α)",
       "log₂ μ = β₀ + β₁·condition + β₂·batch", "Var = μ + α·μ²",
       [(0.95, "β₀", "baseline", False), (2.45, "β₁", "= log₂ fold change", True),
        (3.95, "β₂", "batch effect", False)],
       [(6.45, "μ", "sampling noise\n(Poisson)"), (8.65, "α·μ²", "biological variation\n(α per gene)")],
       "left to the NB:  μ + α·μ²",
       note="The LFC (β₁) comes from the mean model.  How far to trust it (its SE, and so the p-value)\n"
            "comes from the variance model, so DESeq2 has to estimate both.")
print("ok h30 model parts")
