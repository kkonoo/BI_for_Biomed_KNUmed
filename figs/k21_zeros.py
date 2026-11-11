from gallery import *
rng = np.random.default_rng(21)

# Why most of the matrix is zero: only ~10-30% of molecules are captured, so low-expressed genes often read 0
fig = plt.figure(figsize=(9.8, 4.3))
gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 0.95], wspace=0.3, left=0.02, right=0.99, top=0.8, bottom=0.14)
DARK = "#5E6670"; GREYF = "#C9CED6"
# (a) capture
a = fig.add_subplot(gs[0]); a.set_xlim(0, 10); a.set_ylim(0, 10); a.axis("off")
a.set_title("(a) Capture is a lottery", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
a.text(0, 1.03, "each molecule caught with p ≈ 0.2", transform=a.transAxes, fontsize=9, color=TEXT2)
for k, (g, m, y) in enumerate((("gene A", 20, 7.4), ("gene B", 2, 3.6))):
    caught = np.zeros(m, bool); caught[rng.choice(m, 4 if m == 20 else 0, replace=False)] = True
    a.text(0, y + 1.25, f"{g}: {m} molecules in the cell", ha="left", va="center", fontsize=9, color=TEXT, fontweight="semibold")
    for i in range(m):
        a.add_patch(mp.Circle((0.35 + (i % 10) * 0.75, y - (i // 10) * 0.75), 0.24, fc=CRIMSON if caught[i] else "white",
                              ec=CRIMSON if caught[i] else GREYF, lw=1))
    a.text(0, y - 1.5 if m == 20 else y - 0.9, f"counted: {caught.sum()}", ha="left", va="center", fontsize=9.5,
           color=CRIMSON if caught.sum() else TEXT, fontweight="bold")
a.text(0, 1.6, "gene B is expressed, but the matrix says 0", ha="left", va="center", fontsize=8.6, color=TEXT2)
# (b) P(zero | m molecules)
b = fig.add_subplot(gs[1]); m = np.arange(0, 31)
for p, c, lab in ((0.1, CRIMSON, "capture 10%"), (0.3, DARK, "capture 30%")):
    b.plot(m, (1 - p) ** m, color=c, lw=2)
    b.text(m[12] + 0.6, (1 - p) ** 12 + 0.04, lab, fontsize=8.6, color=c, fontweight="semibold")
b.set_xlim(0, 30); b.set_ylim(0, 1.02)
b.set_xlabel("molecules of the gene in the cell", color=TEXT2, fontsize=9)
b.set_ylabel("chance the count is 0", color=TEXT2, fontsize=9)
b.set_yticks([0, 0.5, 1]); b.set_yticklabels(["0", "50%", "100%"])
for sp in ("left", "bottom"): b.spines[sp].set_color(LINE)
b.tick_params(colors=TEXT2, labelsize=8.3)
b.set_title("(b) Low expression → mostly zeros", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
b.text(0, 1.03, "a zero alone says almost nothing", transform=b.transAxes, fontsize=9, color=TEXT2)
# (c) a slice of a real-looking matrix
c = fig.add_subplot(gs[2])
mu = np.exp(rng.normal(-2.2, 1.4, 40)); sf = np.exp(rng.normal(0, 0.4, 30))
M = rng.poisson(np.outer(sf, mu) * 0.35)
c.imshow(np.where(M > 0, np.log1p(M), np.nan), cmap=plt.matplotlib.colors.LinearSegmentedColormap.from_list(
    "z", [tint(0.3), CRIMSON]), aspect="auto", interpolation="nearest")
c.set_facecolor("#F4F5F7"); c.set_xticks([]); c.set_yticks([])
for sp in c.spines.values(): sp.set_visible(False)
c.set_xlabel("40 genes", color=TEXT2, fontsize=9); c.set_ylabel("30 cells", color=TEXT2, fontsize=9)
c.set_title(f"(c) {100 * (M == 0).mean():.0f}% of entries are 0", loc="left", fontsize=12.3, fontweight="semibold",
            color=TEXT, pad=22)
c.text(0, 1.03, "grey = 0, red = detected", transform=c.transAxes, fontsize=9, color=TEXT2)
header(a, "Single-cell RNA-seq · sparsity", "Most zeros mean 'not caught', not 'not expressed'", dy=40, x1=3.55)
fig.savefig(OUT + "sparsity.png"); print("ok k21 zeros", (M == 0).mean())
