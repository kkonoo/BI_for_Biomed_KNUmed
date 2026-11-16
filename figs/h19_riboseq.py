from gallery import *

# Ribosome profiling: freeze ribosomes, digest unprotected RNA, sequence ~30 nt footprints → translation per gene
W, H = 9.8, 4.3
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, (H - 1.0) / H]); ax.set_xlim(0, W); ax.set_ylim(0, H - 1.0); ax.axis("off")
BG = "#F8F9FA"; ARROW = "#9AA1AB"; DARK = "#5E6670"; GREYF = "#C3C8D0"
steps = [("Freeze ribosomes", "a drug stops them where they are"),
         ("Digest bare RNA", "RNase cuts what no ribosome covers"),
         ("Count footprints", "per gene, next to RNA-seq")]
CW, GAP, TOP, CH = 3.08, 0.2, H - 1.12, 2.45
X0 = [0.08 + k * (CW + GAP) for k in range(3)]
for k, (name, sub) in enumerate(steps):
    x = X0[k]
    ax.add_patch(mp.FancyBboxPatch((x, TOP - CH), CW, CH, boxstyle="round,pad=0,rounding_size=0.12",
                                   fc=BG, ec=LINE, lw=0.9))
    ax.scatter([x + 0.27], [TOP - 0.28], s=300, color=tint(0.16), lw=0, zorder=3)
    ax.text(x + 0.27, TOP - 0.28, str(k + 1), ha="center", va="center", fontsize=9.5, fontweight="bold",
            color=CRIMSON, zorder=4)
    ax.text(x + 0.5, TOP - 0.28, name, ha="left", va="center", fontsize=11, fontweight="bold", color=TEXT)
    ax.text(x + 0.16, TOP - 0.6, sub, ha="left", va="center", fontsize=8.6, color=TEXT2)
    if k < 2:
        ax.annotate("", (X0[k + 1] - 0.02, TOP - CH / 2), (x + CW + 0.02, TOP - CH / 2),
                    arrowprops=dict(arrowstyle="-|>", color=ARROW, lw=1.1, mutation_scale=9))

def ribosome(cx, y, a=1.0):
    ax.add_patch(mp.Ellipse((cx, y + 0.11), 0.36, 0.22, fc=DARK, ec="white", lw=0.8, alpha=a, zorder=4))
    ax.add_patch(mp.Ellipse((cx, y + 0.3), 0.26, 0.17, fc=DARK, ec="white", lw=0.8, alpha=a, zorder=4))

RIBO = (1.08, 1.53, 1.98, 2.43)
# 1 — mRNA with ribosomes on the coding region
x = X0[0]; my = TOP - 1.45
ax.plot([x + 0.25, x + 2.85], [my, my], color=DARK, lw=1.4)
ax.add_patch(mp.Rectangle((x + 0.85, my - 0.06), 1.8, 0.12, fc=DARK, ec="none"))
for r in RIBO: ribosome(x + r, my)
ax.text(x + 0.25, my - 0.25, "5′", ha="left", va="center", fontsize=8, color=TEXT2)
ax.text(x + 2.85, my - 0.25, "3′", ha="right", va="center", fontsize=8, color=TEXT2)
ax.text(x + 1.75, my - 0.25, "coding region", ha="center", va="center", fontsize=8, color=TEXT2)
ax.text(x + 0.2, TOP - CH + 0.2, "e.g. cycloheximide", ha="left", va="center", fontsize=8.3, color=TEXT2)

# 2 — only the ~30 nt under each ribosome survives
x = X0[1]
ax.plot([x + 0.25, x + 2.85], [my, my], color=GREYF, lw=1.1, ls=(0, (2, 2)))
for r in RIBO:
    ribosome(x + r, my, a=0.35)
    ax.plot([x + r - 0.15, x + r + 0.15], [my, my], color=CRIMSON, lw=4, solid_capstyle="butt", zorder=5)
ax.text(x + 1.75, my - 0.3, "footprints (~30 nt)", ha="center", va="center", fontsize=8.4, color=CRIMSON,
        fontweight="semibold")
ax.text(x + 0.2, TOP - CH + 0.2, "sequence the footprints, not the mRNA", ha="left", va="center", fontsize=8.3,
        color=TEXT)

# 3 — same mRNA, different footprints → translation efficiency
x = X0[2]; base = TOP - 2.0; sc = 0.4
for g, (name, rna, rib) in enumerate((("gene P", 1.8, 2.2), ("gene Q", 1.8, 0.5))):
    gx = x + 0.45 + g * 1.4
    ax.add_patch(mp.Rectangle((gx, base), 0.36, rna * sc, fc=GREYF, ec="none"))
    ax.add_patch(mp.Rectangle((gx + 0.44, base), 0.36, rib * sc, fc=CRIMSON if g == 0 else tint(0.45), ec="none"))
    ax.text(gx + 0.4, base - 0.14, name, ha="center", va="center", fontsize=8.6, fontweight="bold", color=TEXT)
    ax.text(gx + 0.4, base - 0.33, "TE high" if g == 0 else "TE low", ha="center", va="center", fontsize=7.8,
            color=CRIMSON if g == 0 else TEXT2)
ax.plot([x + 0.3, x + 2.95], [base, base], color=LINE, lw=1)
ax.scatter([x + 0.5], [TOP - 0.93], marker="s", s=22, color=GREYF)
ax.text(x + 0.62, TOP - 0.93, "RNA-seq", ha="left", va="center", fontsize=8, color=TEXT2)
ax.scatter([x + 1.45], [TOP - 0.93], marker="s", s=22, color=CRIMSON)
ax.text(x + 1.57, TOP - 0.93, "Ribo-seq (footprints)", ha="left", va="center", fontsize=8, color=CRIMSON)

bx0, bx1, by = 0.08, W - 0.08, 0.12
ax.add_patch(mp.FancyBboxPatch((bx0, by), bx1 - bx0, 0.46, boxstyle="round,pad=0,rounding_size=0.08",
                               fc=CRIMSON, ec="none"))
ax.text((bx0 + bx1) / 2, by + 0.23, "Translation efficiency (TE) = Ribo-seq ÷ RNA-seq: how hard each mRNA is being translated",
        ha="center", va="center", fontsize=10.3, fontweight="bold", color="white")
header(ax, "Beyond RNA-seq · translation (ribosome profiling)", "Sequence only the RNA that ribosomes are sitting on")
fig.savefig(OUT + "ribo_seq.png"); print("ok h19 ribo-seq")
