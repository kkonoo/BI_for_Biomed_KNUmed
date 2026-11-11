from gallery import *

# One bulk 2× result, three biologies: every cell shifts, a few cells shift a lot, or only the mixture changes
W, H = 9.8, 4.55
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, (H - 0.8) / H]); ax.set_xlim(0, W); ax.set_ylim(0, H - 0.8); ax.axis("off")
DARK = "#5E6670"
def col(v): return tint(min(1.0, 0.1 + v / 6 * 0.9))

# bulk banner
ax.add_patch(mp.FancyBboxPatch((0.1, 3.27), W - 0.2, 0.44, boxstyle="round,pad=0,rounding_size=0.06",
                               fc=FAINT, ec="none"))
ax.text(0.3, 3.49, "Bulk RNA-seq:  gene X  2.0× up in disease,  p = 0.001", ha="left", va="center",
        fontsize=11.5, fontweight="bold", color=TEXT)
ax.text(W - 0.3, 3.49, "all three panels below give exactly this", ha="right", va="center", fontsize=9,
        color=TEXT2)

L, Hh = "o", "s"
cases = [("① Every cell shifts a little", [(1.0, L)] * 10, [(2.0, L)] * 10,
          "every cell: 1 → 2", "a modest, general response", False),
         ("② A few cells shift a lot", [(1.0, L)] * 10, [(1.0, L)] * 8 + [(6.0, L)] * 2,
          "2 of 10 cells: 1 → 6", "usually the interesting case;\ncan be diluted away in bulk", True),
         ("③ Only the mixture changed", [(3.0, Hh)] * 2 + [(0.5, L)] * 8, [(3.0, Hh)] * 6 + [(0.5, L)] * 4,
          "no cell changed; type H 20% → 60%", "a composition shift,\nnot regulation", True)]
PW = 3.1
for k, (title, ctrl, dis, what, means, hl) in enumerate(cases):
    x0 = 0.15 + k * (PW + 0.17)
    ax.text(x0, 2.95, title, ha="left", va="center", fontsize=11.3, fontweight="semibold",
            color=CRIMSON if k == 1 else TEXT)
    for r, (lab, cells) in enumerate((("control", ctrl), ("disease", dis))):
        y = 2.37 - r * 0.55
        ax.text(x0, y + 0.25, lab, ha="left", va="center", fontsize=8.3, color=TEXT2)
        for i, (v, m) in enumerate(cells):
            ax.scatter([x0 + 0.12 + i * 0.235], [y - 0.02], s=150 if m == L else 120, marker=m, color=col(v),
                       lw=0.6, edgecolor="white" if m == L else DARK, zorder=3)
        mean = np.mean([v for v, _ in cells])
        ax.text(x0 + 2.55, y - 0.02, f"mean {mean:.1f}", ha="left", va="center", fontsize=8.6,
                color=TEXT if r else TEXT2, fontweight="bold" if r else "normal")
    ax.plot([x0, x0 + PW - 0.05], [1.53, 1.53], color=LINE, lw=0.8)
    ax.text(x0, 1.33, what, ha="left", va="center", fontsize=9.3, color=TEXT, fontweight="semibold")
    ax.text(x0, 1.13, means, ha="left", va="top", fontsize=9, color=CRIMSON if hl else TEXT2,
            linespacing=1.25)
# legend
ax.text(0.15, 0.55, "colour = expression of gene X in that cell", ha="left", va="center", fontsize=8.5,
        color=TEXT2)
for i, v in enumerate(np.linspace(0, 6, 7)):
    ax.scatter([3.0 + i * 0.16], [0.55], s=70, color=col(v), lw=0, zorder=3)
ax.text(2.88, 0.55, "low", ha="right", va="center", fontsize=8, color=TEXT2)
ax.text(4.08, 0.55, "high", ha="left", va="center", fontsize=8, color=TEXT2)
ax.scatter([6.62], [0.55], s=60, marker="s", color=tint(0.5), edgecolor=DARK, lw=0.6)
ax.text(6.74, 0.55, "cell type H (expresses X highly)", ha="left", va="center", fontsize=8.5, color=TEXT2)
header(ax, "Single-cell RNA-seq · why single cell", "One bulk result, three different biologies", dy=0)
fig.savefig(OUT + "bulk_vs_singlecell.png"); print("ok k1 bulk")
