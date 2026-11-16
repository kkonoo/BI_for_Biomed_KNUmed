from gallery import *

# Nascent transcription (GRO-seq, TT-seq): label RNA made in the last minutes, pull it out, read synthesis
W, H = 9.8, 4.3
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, (H - 1.0) / H]); ax.set_xlim(0, W); ax.set_ylim(0, H - 1.0); ax.axis("off")
BG = "#F8F9FA"; ARROW = "#9AA1AB"; DARK = "#5E6670"; GREYF = "#C3C8D0"
steps = [("Label new RNA", "a pulse of minutes"), ("Pull out the label", "keep only new RNA"),
         ("Read synthesis", "not the accumulated level")]
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

# 1 — Pol II on the gene; the stretch made during the pulse carries the label
x = X0[0]; gy = TOP - 1.55
ax.plot([x + 0.2, x + 2.9], [gy, gy], color=DARK, lw=1.4)
ax.plot([x + 0.2, x + 2.9], [gy - 0.07, gy - 0.07], color=DARK, lw=1.4)
ax.text(x + 0.2, gy - 0.25, "gene (DNA)", ha="left", va="center", fontsize=8, color=TEXT2)
px = x + 2.1
ax.add_patch(mp.Ellipse((px, gy - 0.03), 0.5, 0.36, fc=DARK, ec="none", zorder=4))
ax.text(px, gy - 0.03, "Pol II", ha="center", va="center", fontsize=7.4, color="white", fontweight="bold", zorder=5)
t = np.linspace(0, 1, 50)
xs = px - 0.1 - 1.55 * t; ys = gy + 0.2 + 0.55 * t ** 1.4
cut = 18
ax.plot(xs[:cut], ys[:cut], color=CRIMSON, lw=3, solid_capstyle="butt", zorder=3)
ax.plot(xs[cut - 1:], ys[cut - 1:], color=GREYF, lw=3, solid_capstyle="butt", zorder=3)
ax.text(px - 0.05, gy + 0.32, "made during\nthe pulse", ha="left", va="bottom", fontsize=8, color=CRIMSON,
        fontweight="semibold", linespacing=1.2)
ax.text(x + 0.2, TOP - CH + 0.42, "GRO-seq: BrU, in isolated nuclei", ha="left", va="center", fontsize=8.3,
        color=TEXT)
ax.text(x + 0.2, TOP - CH + 0.2, "TT-seq: 4sU, in living cells", ha="left", va="center", fontsize=8.3, color=TEXT)

# 2 — bead keeps labelled RNA; old RNA is washed away
x = X0[1]; cy = TOP - 1.3
ax.add_patch(mp.Circle((x + 2.35, cy), 0.28, fc=tint(0.2), ec=CRIMSON, lw=1.1, zorder=3))
ax.text(x + 2.35, cy, "bead", ha="center", va="center", fontsize=7.4, color=CRIMSON, zorder=4)
for dy in (0.17, -0.02, -0.2):
    ax.plot([x + 1.35, x + 2.07], [cy + dy, cy + dy * 0.4], color=CRIMSON, lw=3, solid_capstyle="butt", zorder=3)
for i, (xx, yy) in enumerate(((0.25, -0.5), (0.75, -0.66), (1.3, -0.55))):
    ax.plot([x + xx, x + xx + 0.6], [cy + yy, cy + yy], color=GREYF, lw=3, solid_capstyle="butt")
ax.text(x + 0.25, cy + 0.25, "new RNA\n(labelled)", ha="left", va="center", fontsize=8, color=CRIMSON,
        fontweight="semibold", linespacing=1.2)
ax.text(x + 2.1, cy - 0.6, "old RNA\nwashed away", ha="left", va="center", fontsize=8, color=TEXT2,
        linespacing=1.2)
ax.text(x + 0.2, TOP - CH + 0.2, "sequence only what was pulled out", ha="left", va="center", fontsize=8.3,
        color=TEXT)

# 3 — gene A (transcription ×2) and gene B (half-life ×2) from the water-level figure
x = X0[2]; base = TOP - 2.0; sc = 0.42
groups = (("gene A", "transcription ×2", 2.0, 2.0), ("gene B", "half-life ×2", 2.0, 1.0))
for g, (name, why, rs, ns) in enumerate(groups):
    gx = x + 0.45 + g * 1.4
    ax.add_patch(mp.Rectangle((gx, base), 0.36, rs * sc, fc=GREYF, ec="none"))
    ax.add_patch(mp.Rectangle((gx + 0.44, base), 0.36, ns * sc, fc=CRIMSON if g == 0 else tint(0.45), ec="none"))
    ax.text(gx + 0.4, base - 0.14, name, ha="center", va="center", fontsize=8.6, fontweight="bold", color=TEXT)
    ax.text(gx + 0.4, base - 0.33, why, ha="center", va="center", fontsize=7.8, color=TEXT2)
ax.plot([x + 0.3, x + 2.95], [base, base], color=LINE, lw=1)
ax.scatter([x + 0.5], [TOP - 0.93], marker="s", s=22, color=GREYF); ax.text(x + 0.62, TOP - 0.93, "RNA-seq: same",
                                                                           ha="left", va="center", fontsize=8, color=TEXT2)
ax.scatter([x + 1.75], [TOP - 0.93], marker="s", s=22, color=CRIMSON); ax.text(x + 1.87, TOP - 0.93,
                                                                              "nascent: differs", ha="left",
                                                                              va="center", fontsize=8, color=CRIMSON)

bx0, bx1, by = 0.08, W - 0.08, 0.12
ax.add_patch(mp.FancyBboxPatch((bx0, by), bx1 - bx0, 0.46, boxstyle="round,pad=0,rounding_size=0.08",
                               fc=CRIMSON, ec="none"))
ax.text((bx0 + bx1) / 2, by + 0.23, "Nascent RNA measures the inflow (synthesis), so A and B are no longer the same",
        ha="center", va="center", fontsize=10.3, fontweight="bold", color="white")
header(ax, "Beyond RNA-seq · nascent transcription (GRO-seq, TT-seq)", "Catch RNA while it is being made")
fig.savefig(OUT + "nascent_rna.png"); print("ok h17 nascent")
