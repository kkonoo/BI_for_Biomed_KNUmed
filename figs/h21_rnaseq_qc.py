from gallery import *

# RNA-seq QC: where the reads went (per sample), and gene-body coverage (intact vs degraded RNA)
fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.2), gridspec_kw={"wspace": 0.28, "width_ratios": [1.15, 1]})
GREYF = "#C3C8D0"; DARK = "#5E6670"

def style(ax, title, sub):
    for sp in ("left", "bottom"): ax.spines[sp].set_color(LINE)
    ax.tick_params(colors=TEXT2, labelsize=8.5)
    ax.set_title(title, loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=24)
    ax.text(0, 1.03, sub, transform=ax.transAxes, ha="left", va="bottom", fontsize=9.1, color=TEXT2)

# (a) read fates per sample — S5 is the bad one
ax = axes[0]
samples = ["S1", "S2", "S3", "S4", "S5", "S6"]
parts = [("exonic", tint(0.75)), ("intronic", tint(0.35)), ("intergenic", tint(0.15)),
         ("rRNA", DARK), ("unmapped", GREYF)]
vals = np.array([[72, 12, 4, 2, 10], [70, 13, 4, 3, 10], [74, 11, 3, 2, 10],
                 [71, 12, 4, 3, 10], [38, 9, 4, 27, 22], [73, 11, 4, 2, 10]], float)
for i, s in enumerate(samples):
    left = 0
    for j, (lab, col) in enumerate(parts):
        ax.barh(i, vals[i, j], left=left, color=col, height=0.62, lw=0)
        left += vals[i, j]
ax.set_yticks(range(6)); ax.set_yticklabels(samples)
for t in ax.get_yticklabels():
    if t.get_text() == "S5": t.set_color(CRIMSON); t.set_fontweight("bold")
ax.invert_yaxis(); ax.set_xlim(0, 138); ax.set_xticks([0, 50, 100]); ax.spines["bottom"].set_bounds(0, 100)
ax.set_xlabel("% of reads", color=TEXT2, fontsize=9, x=100 / 138 / 2)
ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)
x = 0
for lab, col in parts:
    ax.scatter([x + 2], [-0.95], marker="s", s=24, color=col, clip_on=False)
    ax.text(x + 5, -0.95, lab, ha="left", va="center", fontsize=8.3, color=TEXT2)
    x += 9 + len(lab) * 2.45
ax.text(103, 4, "high rRNA,\nlow mapping\n→ check or drop", ha="left", va="center",
        fontsize=8.6, color=CRIMSON, fontweight="semibold", linespacing=1.3)
style(ax, "(a) Where did the reads go?", "one bar per sample (MultiQC, Qualimap, RSeQC)")

# (b) gene-body coverage: degraded RNA piles up at the 3′ end in a poly(A) library
ax = axes[1]
p = np.linspace(0, 100, 101)
good = 0.72 + 0.22 * np.sin(np.pi * p / 100) ** 0.6
bad = 0.1 + 0.85 * (p / 100) ** 2.4
ax.plot(p, good, color=DARK, lw=2.2); ax.plot(p, bad, color=CRIMSON, lw=2.2)
ax.text(50, 1.0, "intact RNA (RIN 9): even", ha="center", va="center", fontsize=9, color=DARK, fontweight="semibold")
ax.text(30, 0.22, "degraded RNA (RIN 4)", ha="left", va="center", fontsize=9, color=CRIMSON, fontweight="semibold")
ax.text(97, 0.45, "3′ bias", ha="right", va="center", fontsize=9, color=CRIMSON)
ax.set_xlim(0, 100); ax.set_ylim(0, 1.12); ax.set_yticks([])
ax.set_xticks([0, 50, 100]); ax.set_xticklabels(["5′ end", "middle", "3′ end"])
ax.set_xlabel("position along each transcript", color=TEXT2, fontsize=9)
ax.set_ylabel("read coverage", color=TEXT2, fontsize=9)
style(ax, "(b) Gene-body coverage", "averaged over all genes in one sample")

header(axes[0], "Pre-processing · QC", "Check every sample before you trust its counts", dy=44, x1=2.3)
fig.savefig(OUT + "rnaseq_qc.png"); print("ok h21 qc")
