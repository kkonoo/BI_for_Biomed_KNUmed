from gallery import *
rng = np.random.default_rng(9)

# Barcode rank (knee) plot: most barcodes are empty droplets; a hard knee cut loses small real cells
fig = plt.figure(figsize=(9.8, 4.7))
gs = fig.add_gridspec(1, 2, width_ratios=[1.45, 1], wspace=0.22, left=0.07, right=0.99, top=0.8, bottom=0.13)
a = fig.add_subplot(gs[0]); DARK = "#5E6670"; GREYF = "#B8BEC7"
big = np.exp(rng.normal(np.log(6000), 0.55, 6000))
small = np.exp(rng.normal(np.log(420), 0.35, 1200))
empty = np.concatenate([np.exp(rng.normal(np.log(30), 0.75, 150000)), np.exp(rng.normal(np.log(260), 0.4, 2500))])
umi = np.concatenate([big, small, empty]); cell = np.r_[np.ones(big.size + small.size), np.zeros(empty.size)]
o = np.argsort(-umi); umi, cell = umi[o], cell[o]; rank = np.arange(1, umi.size + 1)
a.plot(rank, umi, color=GREYF, lw=2.2, zorder=2)
ci = np.where(cell == 1)[0]
a.scatter(rank[ci], umi[ci], s=3, color=CRIMSON, lw=0, zorder=3)
knee = 6300
a.axvline(knee, color=TEXT, lw=0.9, ls=(0, (4, 3)))
lost = int(((cell == 1) & (rank > knee)).sum())
a.set_xscale("log"); a.set_yscale("log"); a.set_xlim(1, 3e5); a.set_ylim(3, 1.5e5)
a.text(knee * 1.12, 6e4, "hard cut\nat the knee", ha="left", va="center", fontsize=8.5, color=TEXT, linespacing=1.15)
a.text(3, 4.5e3, "① plateau: real cells", ha="left", va="center", fontsize=9, color=CRIMSON, fontweight="semibold")
a.text(3, 2.4e3, "thousands of UMIs each", ha="left", va="center", fontsize=8.4, color=TEXT2)
a.annotate(f"② the cliff: small, low-RNA cells\n    mixed with RNA-rich empty droplets\n    → {lost:,} real cells right of the cut",
           (9000, 330), xytext=(4, 230), textcoords="data", va="center", fontsize=8.4, color=TEXT, linespacing=1.25,
           arrowprops=dict(arrowstyle="-|>", color=TEXT2, lw=0.8, mutation_scale=8))
a.annotate("③ long tail: empty droplets\n    ambient RNA only, ~10–100 UMIs\n    the vast majority of barcodes",
           (6e4, 22), xytext=(4, 11), textcoords="data", va="center", fontsize=8.4, color=TEXT2, linespacing=1.25,
           arrowprops=dict(arrowstyle="-|>", color=TEXT2, lw=0.8, mutation_scale=8))
a.set_xlabel("barcodes, ranked by total UMI", color=TEXT2, fontsize=9.3)
a.set_ylabel("total UMI per barcode", color=TEXT2, fontsize=9.3)
for sp in ("left", "bottom"): a.spines[sp].set_color(LINE)
a.tick_params(colors=TEXT2, labelsize=8.3); a.minorticks_off()
a.set_title("(a) Barcode rank plot (knee plot)", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
a.text(0, 1.03, "one 10x channel · red = barcodes that really hold a cell", transform=a.transAxes, fontsize=9, color=TEXT2)

# (b) EmptyDrops: profile, not total
b = fig.add_subplot(gs[1]); b.set_xlim(0, 10); b.set_ylim(0.2, 9.6); b.axis("off")
b.set_title("(b) EmptyDrops: profile, not total", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
b.text(0, 1.03, "is this barcode different from the soup?", transform=b.transAxes, fontsize=9, color=TEXT2)
genes = ["HBB", "ALB", "MALAT1", "CD3E", "IL7R"]
soup = np.array([0.45, 0.3, 0.2, 0.03, 0.02])
rows = [("soup profile", "learned from the tail ③", soup, DARK, None),
        ("barcode A · 320 UMI", "", np.array([0.08, 0.05, 0.22, 0.35, 0.3]), CRIMSON, "differs\n→ keep (a small cell)"),
        ("barcode B · 410 UMI", "", np.array([0.42, 0.33, 0.19, 0.04, 0.02]), GREYF, "same as soup\n→ drop")]
for k, (t, sub, v, c, verdict) in enumerate(rows):
    y0 = 7.0 - k * 3.4
    b.text(0, y0 + 2.05, t, ha="left", va="center", fontsize=9.6, fontweight="semibold", color=TEXT)
    if sub: b.text(5.3, y0 + 2.05, sub, ha="left", va="center", fontsize=8.2, color=TEXT2)
    for i, val in enumerate(v):
        b.add_patch(mp.Rectangle((0.1 + i * 1.15, y0), 0.8, val * 3.0, fc=c, ec="none"))
        if k == 0: b.text(0.5 + i * 1.15, y0 - 0.35, genes[i], ha="center", va="center", fontsize=6.8, color=TEXT2)
    b.plot([0, 5.9], [y0, y0], color=LINE, lw=0.8)
    if verdict:
        b.text(6.3, y0 + 0.8, verdict, ha="left", va="center", fontsize=8.6, linespacing=1.2,
               color=CRIMSON if c == CRIMSON else TEXT2, fontweight="semibold" if c == CRIMSON else "normal")
header(a, "Single-cell RNA-seq · empty droplets", "Most barcodes are not cells", dy=40, x1=1.86)
fig.savefig(OUT + "knee_plot.png"); print("ok k9 knee", lost)
