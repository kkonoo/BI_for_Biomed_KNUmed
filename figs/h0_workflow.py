from gallery import *

# The bench workflow of bulk RNA-seq: four steps, and the decision made at each one
W, H = 9.8, 4.6
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, (H - 1.0) / H]); ax.set_xlim(0, W); ax.set_ylim(0, H - 1.0); ax.axis("off")
BG = "#F8F9FA"; ARROW = "#9AA1AB"
steps = [
    ("Sample", "tissue, cells, blood",
     ["how many biological replicates?", "spread groups over batches", "record age, sex, time"],
     "§2.4–2.5"),
    ("RNA extraction", "RNA degrades within minutes",
     ["work fast and cold", "measure quality (RIN)", "same kit, same person if possible"],
     "§2.5"),
    ("Library prep", "RNA → cDNA with adapters",
     ["poly(A) or rRNA depletion", "stranded or not", "RNA input, kit"],
     "§2.1–2.2"),
    ("Sequencing", "reads per sample",
     ["read length, SE or PE", "depth (often 20–30 M)", "spread groups over runs"],
     "§2.3"),
]
CW, GAP, TOP, CH = 2.18, 0.28, H - 1.15, 2.45
for k, (name, sub, items, ref) in enumerate(steps):
    x = 0.08 + k * (CW + GAP)
    ax.add_patch(mp.FancyBboxPatch((x, TOP - CH), CW, CH, boxstyle="round,pad=0,rounding_size=0.12",
                                   fc=BG, ec=LINE, lw=0.9))
    ax.scatter([x + 0.3], [TOP - 0.3], s=330, color=tint(0.16), lw=0, zorder=3)
    ax.text(x + 0.3, TOP - 0.3, str(k + 1), ha="center", va="center", fontsize=10, fontweight="bold",
            color=CRIMSON, zorder=4)
    ax.text(x + 0.55, TOP - 0.3, name, ha="left", va="center", fontsize=11.5, fontweight="bold", color=TEXT)
    ax.text(x + 0.18, TOP - 0.62, sub, ha="left", va="center", fontsize=8.8, color=TEXT2)
    ax.text(x + 0.18, TOP - 0.98, "decide here", ha="left", va="center", fontsize=8.6, color=CRIMSON,
            fontweight="semibold")
    for i, it in enumerate(items):
        y = TOP - 1.3 - i * 0.33
        ax.scatter([x + 0.24], [y], marker="s", s=16, color=tint(0.55), lw=0)
        ax.text(x + 0.36, y, it, ha="left", va="center", fontsize=8.6, color=TEXT)
    ax.text(x + CW - 0.15, TOP - CH + 0.17, ref, ha="right", va="center", fontsize=8.4, color=TEXT2)
    if k < 3:
        ax.annotate("", (x + CW + GAP - 0.04, TOP - CH / 2), (x + CW + 0.04, TOP - CH / 2),
                    arrowprops=dict(arrowstyle="-|>", color=ARROW, lw=1.2, mutation_scale=11))
# outcome bar
bx0, bx1, by = 0.08, W - 0.08, 0.12
ax.add_patch(mp.FancyBboxPatch((bx0, by), bx1 - bx0, 0.5, boxstyle="round,pad=0,rounding_size=0.08",
                               fc=CRIMSON, ec="none"))
ax.text((bx0 + bx1) / 2, by + 0.25, "Each choice sets a limit on the analysis — none of them can be undone later",
        ha="center", va="center", fontsize=10.5, fontweight="bold", color="white")
header(ax, "Bulk RNA-seq · the bench workflow", "The analysis starts at the bench", dy=0, x1=1)
fig.savefig(OUT + "rnaseq_workflow.png"); print("ok h0 workflow")
