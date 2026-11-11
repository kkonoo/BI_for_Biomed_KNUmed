from gallery import *

# From association to causation: four explanations for "patients have less Faecalibacterium", checked in order
W, H = 9.8, 4.6
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, (H - 0.8) / H]); ax.set_xlim(0, W); ax.set_ylim(0, H - 0.8); ax.axis("off")
DARK = "#5E6670"; GREYF = "#B8BEC7"
steps = [("1  Artifact", "cases were extracted,\nstored or sequenced\ndifferently",
          "same protocol and batch\nfor both groups; blanks;\nequal depth", True),
         ("2  Confounding", "cases eat less fibre\nor take different\ndrugs",
          "record diet and drugs;\nadjust or match; study\nuntreated patients", True),
         ("3  Reverse causation", "inflammation adds oxygen;\nstrict anaerobes such as\nFaecalibacterium decline",
          "sample people before\nthe disease starts", True),
         ("4  Causation", "less Faecalibacterium\nhelps cause the disease",
          "transfer to germ-free\nmice; intervention\ntrials", False)]
N = len(steps); CW, GAP = 2.18, 0.17
X0 = (W - (N * CW + (N - 1) * GAP)) / 2
for i, (nm, story, check, rule_out) in enumerate(steps):
    x = X0 + i * (CW + GAP)
    col = CRIMSON if rule_out else DARK
    ax.add_patch(mp.FancyBboxPatch((x, 2.55), CW, 0.55, boxstyle="round,pad=0,rounding_size=0.08", fc=col, ec="none"))
    ax.text(x + 0.12, 2.825, nm, ha="left", va="center", fontsize=10, fontweight="bold", color="white")
    if i < N - 1:
        ax.annotate("", (x + CW + GAP, 2.825), (x + CW, 2.825),
                    arrowprops=dict(arrowstyle="-|>", color=GREYF, lw=1, mutation_scale=7, shrinkA=0, shrinkB=0))
    ax.text(x + 0.08, 2.3, "the story", ha="left", va="center", fontsize=7.8, color=TEXT2)
    ax.text(x + 0.08, 2.12, story, ha="left", va="top", fontsize=8.7, color=TEXT, linespacing=1.15)
    ax.plot([x + 0.08, x + CW - 0.08], [1.3, 1.3], color=LINE, lw=0.8)
    ax.text(x + 0.08, 1.12, "how to check" if rule_out else "what would show it", ha="left", va="center", fontsize=7.8,
            color=TEXT2)
    ax.text(x + 0.08, 0.94, check, ha="left", va="top", fontsize=8.7, color=col, linespacing=1.15,
            fontweight="semibold")
xb = X0 + 3 * (CW + GAP) - GAP
ax.plot([X0, xb], [3.22, 3.22], color=CRIMSON, lw=1.4)
ax.text(X0, 3.3, "rule these out first, in this order", ha="left", va="bottom", fontsize=9, color=CRIMSON,
        fontweight="semibold")
ax.text(X0, 0.12, "Patients have less Faecalibacterium than controls. For most such findings, the evidence stops "
        "before step 4.", ha="left", va="center", fontsize=8.8, color=TEXT2)
header(ax, "Metagenomics · causation", "Before calling a difference a cause, rule out three other stories",
       x=X0 / W, dy=0, x1=1 - X0 / W)
fig.savefig(OUT + "causal_order.png"); print("ok n16")
