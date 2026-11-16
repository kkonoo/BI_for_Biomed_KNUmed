from gallery import *

# A DEG is a difference between groups that is large compared with the spread within groups
fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.9), sharey=True, gridspec_kw={"wspace": 0.14})
GREYF = "#9AA1AB"
genes = (("(a) Gene X", "same difference, small spread", [5.0, 5.25, 4.85, 5.1], [7.05, 6.85, 7.2, 6.95],
          "difference ≫ spread → a DEG", CRIMSON),
         ("(b) Gene Y", "same difference, large spread", [3.4, 6.6, 4.4, 5.6], [5.1, 8.6, 6.1, 8.2],
          "difference ≈ spread → cannot tell", TEXT2))
for ax, (title, sub, ctrl, case, msg, mc) in zip(axes, genes):
    for sp in ("left", "bottom"): ax.spines[sp].set_color(LINE)
    ax.tick_params(colors=TEXT2, labelsize=8.5)
    for x0, vals, col in ((0, ctrl, GREYF), (1, case, CRIMSON)):
        xs = x0 + np.array([-0.07, -0.025, 0.025, 0.07])
        ax.scatter(xs, vals, s=55, color=col, edgecolor="white", lw=0.8, zorder=3)
        m = np.mean(vals)
        ax.plot([x0 - 0.18, x0 + 0.18], [m, m], color=TEXT, lw=1.6, zorder=2)
    ax.annotate("", (1.32, np.mean(case)), (1.32, np.mean(ctrl)),
                arrowprops=dict(arrowstyle="<->", color=TEXT, lw=1))
    ax.text(1.36, (np.mean(case) + np.mean(ctrl)) / 2, "mean\ndifference", ha="left", va="center", fontsize=8.4,
            color=TEXT, linespacing=1.2)
    ax.set_xlim(-0.45, 1.75); ax.set_ylim(2.5, 9.6)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["control", "case"])
    ax.text(0.5, 9.2, msg, ha="center", va="center", fontsize=9.2, color=mc, fontweight="semibold")
    ax.set_title(title, loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=24)
    ax.text(0, 1.03, sub, transform=ax.transAxes, ha="left", va="bottom", fontsize=9.1, color=TEXT2)
axes[0].set_ylabel("expression (log scale)", color=TEXT2, fontsize=9)

header(axes[0], "Differential expression · the question", "A DEG is a difference that stands out from the noise",
       dy=44, x1=2.14)
fig.savefig(OUT + "deg_signal_noise.png"); print("ok h24 deg")
