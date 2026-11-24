from gallery import *
rng = np.random.default_rng(9)

# Barcode rank (knee) plot: most barcodes are empty droplets. Cell Ranger calls cells in two steps:
# ① everything above m/10 (m = 99th percentile of the top N barcodes), ② EmptyDrops on the barcodes just below
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
# how Cell Ranger sets the cut: ① m = 99th percentile of the top N barcodes (N = expected cells), cut = m / 10
N = 7000
m = np.percentile(umi[:N], 99); cut = m / 10
r_cut = int((umi >= cut).sum())                     # rank where the curve crosses the cut
r_test = int((umi >= 100).sum())                    # ② below the cut: candidates tested one by one
a.axhline(cut, color=TEXT, lw=0.9, ls=(0, (4, 3)), xmax=0.86)
a.axvspan(r_cut, r_test, color=FAINT, zorder=0)
lost = int(((cell == 1) & (rank > r_cut)).sum())
a.set_xscale("log"); a.set_yscale("log"); a.set_xlim(1, 3e5); a.set_ylim(3, 1.5e5)
a.text(1.3, cut * 1.25, f"① cut = m ÷ 10 = {cut:,.0f} UMI", ha="left", va="bottom",
       fontsize=8.6, color=TEXT, fontweight="semibold")
a.text(1.3, cut * 0.8, "m = UMI of a big cell (99th percentile of the top N expected)\nabove the cut → a cell, no test needed", ha="left", va="top", fontsize=8.2, color=TEXT2, linespacing=1.25)
a.text(3, 1.15e4, "plateau: real cells", ha="left", va="center", fontsize=9, color=CRIMSON, fontweight="semibold")
a.text(3, 7.2e3, "thousands of UMIs each", ha="left", va="center", fontsize=8.4, color=TEXT2)
a.annotate(f"② the knee: small, low-RNA cells mixed\n    with RNA-rich empty droplets\n    → each barcode tested (EmptyDrops, b)\n    {lost:,} real cells sit below the cut",
           (np.sqrt(r_cut * r_test), 300), xytext=(2.2, 95), textcoords="data", va="center", fontsize=8.4, color=TEXT,
           linespacing=1.25, arrowprops=dict(arrowstyle="-|>", color=TEXT2, lw=0.8, mutation_scale=8))
a.annotate("③ long tail: empty droplets, ~10–100 UMIs\n    → never cells; used to learn the soup",
           (6e4, 22), xytext=(2.2, 9), textcoords="data", va="center", fontsize=8.4, color=TEXT2, linespacing=1.25,
           arrowprops=dict(arrowstyle="-|>", color=TEXT2, lw=0.8, mutation_scale=8))
a.set_xlabel("barcodes, ranked by total UMI", color=TEXT2, fontsize=9.3)
a.set_ylabel("total UMI per barcode", color=TEXT2, fontsize=9.3)
for sp in ("left", "bottom"): a.spines[sp].set_color(LINE)
a.tick_params(colors=TEXT2, labelsize=8.3); a.minorticks_off()
a.set_title("(a) Barcode rank plot (knee plot)", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
a.text(0, 1.03, "one 10x channel (simulated) · red = barcodes that really hold a cell", transform=a.transAxes, fontsize=9, color=TEXT2)

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
fig.savefig(OUT + "knee_plot.png"); print("ok k9 knee", lost, round(cut), r_cut, r_test)
