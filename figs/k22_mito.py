from gallery import *

# Why a high mito % flags a broken cell: cytoplasmic mRNA leaks out, mitochondrial transcripts stay protected
W, H = 9.8, 4.2
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, (H - 0.8) / H]); ax.set_xlim(0, W); ax.set_ylim(0, H - 0.8); ax.axis("off")
DARK = "#5E6670"; GREYF = "#B8BEC7"
rng = np.random.default_rng(22)
def cell(cx, cy, broken, n_cyto):
    ls = (0, (3, 2.5)) if broken else "-"
    ax.add_patch(mp.Ellipse((cx, cy), 2.3, 1.85, fc="white", ec=TEXT2, lw=1.4, ls=ls))
    ax.add_patch(mp.Circle((cx - 0.15, cy + 0.05), 0.32, fc=FAINT, ec=LINE, lw=0.8))
    for mx, my in ((cx + 0.6, cy + 0.45), (cx + 0.5, cy - 0.5), (cx - 0.7, cy - 0.5)):
        ax.add_patch(mp.Ellipse((mx, my), 0.42, 0.22, fc=tint(0.15), ec=CRIMSON, lw=0.9))
        ax.plot([mx - 0.1, mx + 0.1], [my, my], color=CRIMSON, lw=1.4)
    pts = []
    while len(pts) < n_cyto:
        x, y = rng.uniform(-1.0, 1.0), rng.uniform(-0.8, 0.8)
        if (x / 1.05) ** 2 + (y / 0.85) ** 2 < 0.85 and np.hypot(x + 0.15, y - 0.05) > 0.4 and \
           min(np.hypot(x - 0.6, y - 0.45), np.hypot(x - 0.5, y + 0.5), np.hypot(x + 0.7, y + 0.5)) > 0.3:
            pts.append((x, y))
    for x, y in pts:
        ax.plot([cx + x - 0.07, cx + x + 0.07], [cy + y, cy + y + 0.03], color=DARK, lw=1.2, solid_capstyle="round")
cell(1.4, 2.0, False, 26); cell(5.3, 2.0, True, 3)
for x, y in ((6.75, 2.7), (6.9, 1.4), (3.85, 2.9), (6.65, 2.1), (3.9, 1.2)):
    ax.plot([x - 0.07, x + 0.07], [y, y + 0.03], color=DARK, lw=1.2, solid_capstyle="round")
ax.text(1.4, 3.3, "Healthy cell", ha="center", va="center", fontsize=11, fontweight="semibold", color=TEXT)
ax.text(5.3, 3.3, "Broken membrane", ha="center", va="center", fontsize=11, fontweight="semibold", color=CRIMSON)
ax.text(5.3, 0.82, "cytoplasmic mRNA leaks out;\nmitochondria keep theirs", ha="center", va="center", fontsize=8.6,
        color=TEXT2, linespacing=1.25)
ax.text(1.4, 0.82, "dark dashes = cytoplasmic mRNA\nred = mitochondrial transcripts", ha="center", va="center",
        fontsize=8.6, color=TEXT2, linespacing=1.25)
# bars
x0 = 7.6
ax.text(x0, 3.3, "What the cell reads", ha="left", va="center", fontsize=11, fontweight="semibold", color=TEXT)
for k, (lab, cyto, mito) in enumerate((("healthy", 970, 30), ("broken", 70, 30))):
    y = 2.45 - k * 1.0; tot = cyto + mito
    ax.text(x0, y + 0.32, f"{lab}: {tot:,} UMI", ha="left", va="center", fontsize=8.8, color=TEXT)
    wbar = 2.0 * tot / 1000
    ax.add_patch(mp.Rectangle((x0, y - 0.13), wbar * cyto / tot, 0.26, fc=DARK, ec="none"))
    ax.add_patch(mp.Rectangle((x0 + wbar * cyto / tot, y - 0.13), wbar * mito / tot, 0.26, fc=CRIMSON, ec="none"))
    ax.text(x0 + wbar + 0.08, y, f"mito {100 * mito / tot:.0f}%", ha="left", va="center", fontsize=9.5,
            fontweight="bold", color=CRIMSON if k else TEXT)
ax.text(x0, 0.82, "same mito reads, fewer others\n→ the ratio rises", ha="left", va="center", fontsize=8.6,
        color=TEXT2, linespacing=1.25)
ax.text(0.15, 0.15, "Read it with the gene count: high mito + few genes = broken.  High mito + normal genes = a mito-rich "
        "cell (cardiomyocyte, hepatocyte).", ha="left", va="center", fontsize=8.8, color=TEXT)
header(ax, "Single-cell RNA-seq · dying cells", "A high mito % is a ratio that rises when the rest leaks out", dy=0)
fig.savefig(OUT + "mito_fraction.png"); print("ok k22 mito")
