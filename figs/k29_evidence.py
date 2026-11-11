from gallery import *

# How far each single-cell result can carry a claim, and where it belongs in a paper
W, H = 9.8, 4.6
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, (H - 0.8) / H]); ax.set_xlim(0, W); ax.set_ylim(0, H - 0.8); ax.axis("off")
DARK = "#5E6670"
rows = [("Clusters with validated markers", 4, "robust · descriptive", "Results"),
        ("Donor-level pseudobulk DE", 4, "robust · inferential", "Results"),
        ("Differential abundance (compositional, n = donors)", 3, "robust if n = donors", "Results, with caveats"),
        ("Pseudotime ordering", 2, "hypothesis", "Results, 'consistent with'"),
        ("Regulon / TF activity / gene programs", 2, "hypothesis", "Results, 'consistent with'"),
        ("RNA velocity direction", 1, "hypothesis · contested method", "Supplementary"),
        ("Ligand–receptor predictions", 1, "hypothesis only", "Discussion or grant aim")]
ax.text(0.15, 3.6, "result", ha="left", va="center", fontsize=8.6, color=TEXT2)
ax.text(4.3, 3.6, "strength", ha="left", va="center", fontsize=8.6, color=TEXT2)
ax.text(7.65, 3.6, "where it goes", ha="left", va="center", fontsize=8.6, color=TEXT2)
for i, (nm, s, lab, where) in enumerate(rows):
    y = 3.2 - i * 0.46
    if i == 3: ax.plot([0.15, W - 0.15], [y + 0.23, y + 0.23], color=CRIMSON, lw=0.9, ls=(0, (4, 3)))
    ax.text(0.15, y, nm, ha="left", va="center", fontsize=9.4, color=TEXT, fontweight="semibold" if s >= 3 else "normal")
    for k in range(4):
        ax.add_patch(mp.FancyBboxPatch((4.3 + k * 0.3, y - 0.09), 0.24, 0.18, boxstyle="round,pad=0,rounding_size=0.04",
                                       fc=(DARK if s >= 3 else tint(0.55)) if k < s else "#E9EBEE", ec="none"))
    ax.text(5.6, y, lab, ha="left", va="center", fontsize=8.5, color=TEXT2 if s >= 3 else CRIMSON)
    ax.text(7.65, y, where, ha="left", va="center", fontsize=8.8, color=TEXT)
ax.text(0.15, 0.02, "Below the dashed line: write 'consistent with', and pair the result with a validation experiment.",
        ha="left", va="center", fontsize=8.6, color=CRIMSON)
header(ax, "Single-cell RNA-seq · strength of evidence", "Robust results vs hypotheses from one snapshot", dy=0)
fig.savefig(OUT + "evidence.png"); print("ok k29 evidence")
