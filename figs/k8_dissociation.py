from gallery import *

# Dissociation is an experiment: fragile cells are lost and survivors switch on stress genes
W, H = 9.8, 4.1
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, (H - 0.8) / H]); ax.set_xlim(0, W); ax.set_ylim(0.3, H - 0.5); ax.axis("off")
DARK = "#5E6670"; BIG = tint(0.45); SMALL = "#8A929C"
big = [(0.75, 2.55, "neuron"), (1.75, 2.75, "adipocyte"), (2.75, 2.5, "cardio-\nmyocyte")]
small = [(0.55 + i * 0.47, 1.65 + (0.12 if i % 2 else 0)) for i in range(6)]

def box(x0, w, title, ec=LINE):
    ax.add_patch(mp.FancyBboxPatch((x0, 0.95), w, 2.55, boxstyle="round,pad=0,rounding_size=0.08",
                                   fc="white", ec=ec, lw=0.9))
    ax.text(x0 + 0.2, 3.3, title, ha="left", va="center", fontsize=11, fontweight="semibold", color=TEXT)

# in the tissue
box(0.2, 3.1, "In the tissue")
for x, y, lab in big:
    ax.add_patch(mp.Circle((x, y), 0.33, fc=BIG, ec="none"))
    ax.text(x, y, lab, ha="center", va="center", fontsize=7.2, color="white", linespacing=1.0)
for x, y in small:
    ax.add_patch(mp.Circle((x + 0.1, y), 0.13, fc=SMALL, ec="none"))
ax.text(0.4, 1.18, "large / fragile 3 · small / robust 6", ha="left", va="center", fontsize=8.6, color=TEXT2)

# the step
ax.annotate("", (4.15, 2.25), (3.4, 2.25), arrowprops=dict(arrowstyle="-|>", color=LINE, lw=1.2, mutation_scale=10))
ax.text(3.75, 2.95, "enzymatic\ndissociation", ha="center", va="center", fontsize=9, fontweight="semibold",
        color=TEXT, linespacing=1.15)
ax.text(3.75, 1.65, "37 °C\n30–60 min\ncells alive", ha="center", va="center", fontsize=8.3, color=TEXT2,
        linespacing=1.25)

# what reaches the sequencer
x0 = 4.25
box(x0, 5.35, "What reaches the sequencer", ec=CRIMSON)
for (x, y, lab) in big:
    ax.add_patch(mp.Circle((x + x0 - 0.2, y), 0.33, fc="none", ec=LINE, lw=1, ls=(0, (2, 2))))
ax.text(x0 + 1.55, 2.0, "destroyed — missing from the data", ha="center", va="center", fontsize=8.3, color=TEXT2)
for x, y in small:
    cx = x + x0 - 0.1
    ax.add_patch(mp.Circle((cx, y - 0.1), 0.13, fc=SMALL, ec="none"))
    ax.add_patch(mp.Circle((cx + 0.1, y - 0.0), 0.045, fc=CRIMSON, ec="none"))
ax.text(x0 + 0.2, 1.18, "large 0 · small 6", ha="left", va="center", fontsize=8.6, color=TEXT2)
# side notes
nx = x0 + 3.45
ax.text(nx, 2.85, "① composition", ha="left", va="center", fontsize=9.5, fontweight="bold", color=CRIMSON)
ax.text(nx, 2.55, "fragile types under-counted", ha="left", va="center", fontsize=8.6, color=TEXT)
ax.text(nx, 2.05, "② the survivors", ha="left", va="center", fontsize=9.5, fontweight="bold", color=CRIMSON)
ax.text(nx, 1.75, "stress genes switched on:", ha="left", va="center", fontsize=8.6, color=TEXT)
ax.text(nx, 1.5, "FOS  JUN  HSPA1A", ha="left", va="center", fontsize=8.6, color=TEXT, fontstyle="italic",
        fontweight="semibold")
ax.text(0.2, 0.55, "Measured composition = composition of the cells that survived.  "
        "Nuclei (snRNA-seq) avoid both problems: they survive freezing and do not respond.",
        ha="left", va="center", fontsize=8.8, color=TEXT2)
header(ax, "Single-cell RNA-seq · cells vs nuclei", "You measure what survived, not what was there", dy=0)
fig.savefig(OUT + "dissociation_bias.png"); print("ok k8 dissociation")
