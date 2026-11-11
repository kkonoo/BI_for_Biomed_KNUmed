from gallery import *

# Between mRNA and protein: translation rate, protein half-life, and modification, none of which RNA-seq sees
fig, axes = plt.subplots(1, 3, figsize=(9.8, 3.9), gridspec_kw={"wspace": 0.28})
DARK = "#5E6670"; GREYF = "#B8BEC7"

def title(ax, t, sub):
    ax.set_title(t, loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
    ax.text(0, 1.03, sub, transform=ax.transAxes, fontsize=9, color=TEXT2)

# (a) translation rate: same mRNA, different ribosome loading
a = axes[0]; a.set_xlim(0, 10); a.set_ylim(0, 10); a.axis("off"); a.set_aspect("equal", adjustable="datalim")
for k, (nrib, nprot, lab) in enumerate(((1, 1, "few ribosomes"), (6, 6, "many ribosomes"))):
    y = 7.6 - k * 4.0
    a.plot([0.3, 5.6], [y, y], color=DARK, lw=2.4, solid_capstyle="round")
    for i in range(nrib):
        x = 0.9 + i * 0.85
        a.add_patch(mp.Ellipse((x, y + 0.42), 0.72, 0.56, fc=GREYF, ec="white", lw=0.8))
    a.text(0.3, y - 0.8, lab, ha="left", va="center", fontsize=8.8, color=TEXT2)
    a.annotate("", (6.7, y), (5.9, y), arrowprops=dict(arrowstyle="-|>", color=GREYF, lw=1, mutation_scale=8))
    for i in range(nprot):
        a.add_patch(mp.Circle((7.15 + (i % 3) * 0.55, y + 0.3 - (i // 3) * 0.55), 0.21,
                              fc=CRIMSON if k else tint(0.45), ec="none"))
a.text(0.3, 0.6, "the same mRNA count, ≥ 100-fold\ndifference in protein made per mRNA", ha="left", va="center",
       fontsize=8.6, color=TEXT2, linespacing=1.2)
title(a, "(a) Translation rate", "mRNA (line) with ribosomes on it")

# (b) protein half-life: same synthesis, different decay
b = axes[1]
t = np.linspace(0, 96, 300)
ks = 1.0
for hl, col, lab in ((48, CRIMSON, "half-life 48 h"), (4, DARK, "half-life 4 h")):
    kd = np.log(2) / hl
    b.plot(t, ks / kd * (1 - np.exp(-kd * t)), color=col, lw=2.2)
    yend = ks / kd * (1 - np.exp(-kd * 96))
    b.text(94, yend + 3, lab, ha="right", va="bottom", fontsize=8.8, color=col, fontweight="semibold")
b.set_xlim(0, 100); b.set_ylim(0, 70)
b.set_xticks([0, 24, 48, 72, 96]); b.set_yticks([])
b.set_xlabel("hours after synthesis starts", color=TEXT2, fontsize=9)
b.set_ylabel("protein amount", color=TEXT2, fontsize=9)
for sp in ("left", "bottom"): b.spines[sp].set_color(LINE)
b.tick_params(colors=TEXT2, labelsize=8.3)
title(b, "(b) Protein half-life", "same mRNA, same synthesis rate")

# (c) modification: same amount, different activity
c = axes[2]; c.set_xlim(0, 10); c.set_ylim(0, 10); c.axis("off"); c.set_aspect("equal", adjustable="datalim")
for k, (ph, lab) in enumerate(((False, "unmodified"), (True, "phosphorylated"))):
    x = 2.2 + k * 4.6
    for i in range(6):
        cx, cy = x - 1.0 + (i % 3) * 1.0, 6.9 - (i // 3) * 1.05
        c.add_patch(mp.Circle((cx, cy), 0.4, fc=GREYF if not ph else DARK, ec="none"))
        if ph:
            c.add_patch(mp.Circle((cx + 0.36, cy + 0.33), 0.15, fc=CRIMSON, ec="white", lw=0.8))
    c.text(x, 5.0, lab, ha="center", va="center", fontsize=9, color=CRIMSON if ph else TEXT2,
           fontweight="semibold" if ph else "normal")
    c.text(x, 4.2, "inactive" if not ph else "active", ha="center", va="center", fontsize=9.5,
           fontweight="bold", color=CRIMSON if ph else TEXT2)
c.text(0.3, 2.0, "the same amount of protein; the\nactivity is set by the modification\n(also: location, complex partners)",
       ha="left", va="center", fontsize=8.6, color=TEXT2, linespacing=1.2)
title(c, "(c) Modification", "red dot = phosphate group")
header(axes[0], "Proteomics · mRNA to protein", "Three steps after the mRNA that RNA-seq cannot see", dy=40, x1=3.6)
fig.savefig(OUT + "mrna_protein_steps.png"); print("ok m8 steps", round((np.log(2) / 4) / (np.log(2) / 48), 1))
