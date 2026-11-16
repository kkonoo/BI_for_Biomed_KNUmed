from gallery import *

# The design formula turns the sample table into a 0/1 design matrix; each column gets one coefficient
fig = plt.figure(figsize=(9.0, 4.0))
ax = fig.add_axes([0.02, 0.02, 0.96, 0.74]); ax.set_xlim(0, 10); ax.set_ylim(0, 4.6); ax.axis("off")
rows = [("S1", "A", "control"), ("S2", "A", "control"), ("S3", "B", "control"),
        ("S4", "A", "case"), ("S5", "B", "case"), ("S6", "B", "case")]
RH, TOP = 0.46, 3.85

def table(x0, cols, widths, cells, hl=None, head_col=TEXT):
    x = x0
    for c, w in zip(cols, widths):
        ax.text(x + w / 2, TOP + 0.38, c, ha="center", va="center", fontsize=9.2, fontweight="bold", color=head_col)
        x += w
    ax.plot([x0, x0 + sum(widths)], [TOP + 0.16, TOP + 0.16], color=LINE, lw=1)
    for i, r in enumerate(cells):
        x = x0; y = TOP - i * RH - 0.1
        for j, (v, w) in enumerate(zip(r, widths)):
            on = hl is not None and hl(j, v)
            if on:
                ax.add_patch(mp.Rectangle((x + 0.08, y - RH / 2 + 0.05), w - 0.16, RH - 0.1, fc=tint(0.2), ec="none"))
            ax.text(x + w / 2, y, v, ha="center", va="center", fontsize=9.5,
                    color=CRIMSON if on else TEXT, fontweight="bold" if on else "normal")
            x += w

table(0.1, ["sample", "batch", "condition"], [0.85, 0.85, 1.15], rows)
ax.annotate("", (4.75, 2.6), (3.25, 2.6), arrowprops=dict(arrowstyle="-|>", color="#9AA1AB", lw=1.3,
                                                         mutation_scale=12))
ax.text(4.0, 2.95, "design =", ha="center", va="center", fontsize=9, color=TEXT2)
ax.text(4.0, 2.25, "~ batch + condition", ha="center", va="center", fontsize=9.6, family="DejaVu Sans Mono",
        color=TEXT, fontweight="bold")
X = [(r[0], "1", "1" if r[1] == "B" else "0", "1" if r[2] == "case" else "0") for r in rows]
table(4.95, ["", "intercept", "batch B", "condition case"], [0.55, 1.05, 1.05, 1.45], X,
      hl=lambda j, v: j == 3 and v == "1")
for x, lab, col in ((5.5 + 0.525, "β₀", TEXT), (6.55 + 0.525, "β₂", TEXT), (7.6 + 0.725, "β₁", CRIMSON)):
    ax.text(x, TOP - 6 * RH - 0.05, lab, ha="center", va="center", fontsize=11.5, fontweight="bold", color=col)
ax.text(5.0, TOP - 6 * RH - 0.05, "coefficient", ha="left", va="center", fontsize=8.6, color=TEXT2)
ax.text(9.95, 0.6, "control and batch A are the reference (0) · β₁ = log₂FC for case vs control, after batch",
        ha="right", va="center", fontsize=8.8, color=TEXT2)
header(ax, "Differential expression · design matrix", "The design formula turns the sample table into numbers")
fig.savefig(OUT + "design_matrix.png"); print("ok h27 design matrix")
