from gallery import *

# Study design: the donor is the replicate (spend on donors, not cells) and processing day must not equal condition
W, H = 9.8, 5.3
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, (H - 0.8) / H]); ax.set_xlim(0, W); ax.set_ylim(-0.6, 3.9); ax.axis("off")
DARK = "#5E6670"; GREYF = "#B8BEC7"
def person(x, y, c, lab=None, s=1.0):
    ax.add_patch(mp.Circle((x, y + 0.17 * s), 0.09 * s, fc=c, ec="none"))
    ax.add_patch(mp.FancyBboxPatch((x - 0.12 * s, y - 0.16 * s), 0.24 * s, 0.24 * s,
                                   boxstyle="round,pad=0,rounding_size=0.08", fc=c, ec="none"))
    if lab: ax.text(x, y - 0.33 * s, lab, ha="center", va="center", fontsize=7.4, color=TEXT2)
def ptitle(x, t, sub):
    ax.text(x, 3.62, t, ha="left", va="center", fontsize=11.6, fontweight="semibold", color=TEXT)
    ax.text(x, 3.36, sub, ha="left", va="center", fontsize=8.7, color=TEXT2)

# (a) same budget, two designs
ptitle(0.1, "(a) Same budget: 40,000 cells", "n for a patient-vs-control test = donors per group")
designs = [("2 vs 2 donors × 10,000 cells", 2, "n = 2", DARK, 2.4), ("8 vs 8 donors × 2,500 cells", 8, "n = 8", CRIMSON, 1.1)]
for lab, k, n, col, y in designs:
    ax.text(0.1, y + 0.55, lab, ha="left", va="center", fontsize=9.3, color=TEXT, fontweight="semibold")
    for i in range(k): person(0.3 + i * 0.22, y, tint(0.6), s=0.8)
    for i in range(k): person(0.5 + (k + i) * 0.22, y, GREYF, s=0.8)
    ax.text(4.55, y, n, ha="right", va="center", fontsize=13, fontweight="bold", color=col)
ax.text(0.1, 0.1, "red = patients, grey = controls.  More cells per donor\nadd resolution, not power (§6)", ha="left",
        va="center", fontsize=8.5, color=TEXT2, linespacing=1.25)
ax.plot([4.85, 4.85], [-0.3, 3.7], color="#E3E6EB", lw=0.9)

# (b) processing day
x0 = 5.1
ptitle(x0, "(b) Processing day must not equal condition", "batch is fixed at the bench")
def lane(x, y, people, title, verdict, good):
    ax.text(x, y + 0.62, title, ha="left", va="center", fontsize=9.3, fontweight="semibold", color=TEXT)
    for r, row in enumerate(people):
        ax.text(x, y - r * 0.62, f"day {r + 1}", ha="left", va="center", fontsize=8.2, color=TEXT2)
        ax.add_patch(mp.FancyBboxPatch((x + 0.55, y - r * 0.62 - 0.24), 1.6, 0.48, boxstyle="round,pad=0,rounding_size=0.06",
                                       fc=FAINT, ec="none"))
        for i, c in enumerate(row): person(x + 0.85 + i * 0.5, y - r * 0.62 - 0.02, tint(0.6) if c == "P" else GREYF, s=0.85)
    ax.text(x + 2.35, y - 0.31, verdict, ha="left", va="center", fontsize=8.5, linespacing=1.2,
            color=TEXT if good else CRIMSON, fontweight="semibold")
lane(x0, 2.45, [["P", "P", "P"], ["C", "C", "C"]], "Confounded", "✕ day = disease:\n  no analysis can\n  separate them", False)
lane(x0, 0.8, [["P", "C", "P"], ["C", "P", "C"]], "Balanced (or pooled)", "✓ day effect is\n  shared by both\n  groups", True)
ax.text(x0, -0.42, "Pooling several donors per lane, then demultiplexing by genotype, also flags cross-donor doublets",
        ha="left", va="center", fontsize=8.2, color=TEXT2)
header(ax, "Single-cell RNA-seq · study design", "The replicate is the donor, not the cell", dy=0)
fig.savefig(OUT + "donor_design.png"); print("ok k19 donors")
