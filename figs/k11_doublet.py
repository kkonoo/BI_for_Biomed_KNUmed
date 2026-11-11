from gallery import *
from sklearn.neighbors import NearestNeighbors
rng = np.random.default_rng(11)

# Doublet detection: simulate artificial doublets, embed with real cells, score each cell by its share of
# artificial neighbours (pANN), then cut at the expected doublet count
fig = plt.figure(figsize=(9.8, 4.7))
gs = fig.add_gridspec(1, 3, width_ratios=[0.9, 1.15, 1.0], wspace=0.25, left=0.02, right=0.99, top=0.8, bottom=0.14)
DARK = "#5E6670"; TC = "#8A929C"; MC = "#C3C8CF"
# (a) simulate
a = fig.add_subplot(gs[0]); a.set_xlim(0, 10); a.set_ylim(0, 10); a.axis("off")
a.set_title("① Simulate", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
a.text(0, 1.03, "average two real profiles", transform=a.transAxes, fontsize=9, color=TEXT2)
genes = ["CD3E", "IL7R", "LYZ", "CD14"]
prof = [("T cell", np.array([0.9, 0.7, 0, 0]), DARK), ("monocyte", np.array([0, 0, 0.95, 0.8]), MC),
        ("artificial\ndoublet", np.array([0.45, 0.35, 0.47, 0.4]), None)]
for k, (nm, v, c) in enumerate(prof):
    y = 7.2 - k * 3.0
    a.text(0.2, y + 0.6, nm, ha="left", va="center", fontsize=9.2, fontweight="semibold",
           color=CRIMSON if c is None else TEXT, linespacing=1.1)
    for i, val in enumerate(v):
        x = 3.6 + i * 1.5
        if c is None:
            a.add_patch(mp.Rectangle((x, y), 1.0, v[i] * 2, fc="white", ec=CRIMSON, lw=1, ls=(0, (2, 1.5))))
        else:
            a.add_patch(mp.Rectangle((x, y), 1.0, val * 2, fc=c, ec="none"))
        if k == 0: a.text(x + 0.5, 9.6, genes[i], ha="center", va="center", fontsize=7.4, color=TEXT2)
    a.plot([3.4, 9.6], [y, y], color=LINE, lw=0.8)
a.text(6.6, 6.35, "+", ha="center", va="center", fontsize=12, color=TEXT2)
a.text(6.6, 3.55, "average ↓", ha="center", va="center", fontsize=8.5, color=TEXT2)
a.text(0.2, 0.15, "× thousands of random pairs", ha="left", va="center", fontsize=8.5, color=TEXT2)

# (b) embed real + artificial, count neighbours
b = fig.add_subplot(gs[1])
T = rng.normal([-2.2, 0.0], 0.55, (300, 2)); M = rng.normal([2.2, 0.6], 0.55, (200, 2))
D = (T[rng.integers(0, 300, 24)] + M[rng.integers(0, 200, 24)]) / 2 + rng.normal(0, 0.25, (24, 2))
real = np.vstack([T, M, D]); kind = np.r_[np.zeros(300), np.ones(200), 2 * np.ones(24)]
pairs = rng.integers(0, 500, (175, 2))
art = (real[pairs[:, 0]] + real[pairs[:, 1]]) / 2 + rng.normal(0, 0.2, (175, 2))
allp = np.vstack([real, art]); isart = np.r_[np.zeros(len(real)), np.ones(len(art))]
k = 8
nn = NearestNeighbors(n_neighbors=k + 1).fit(allp)
_, idx = nn.kneighbors(real)
pann = isart[idx[:, 1:]].mean(1)
b.scatter(T[:, 0], T[:, 1], s=9, color=DARK, lw=0, alpha=0.8)
b.scatter(M[:, 0], M[:, 1], s=9, color=MC, lw=0, alpha=0.9)
b.scatter(D[:, 0], D[:, 1], s=14, color=CRIMSON, lw=0, zorder=3)
b.scatter(art[:, 0], art[:, 1], s=16, marker="D", facecolor="none", edgecolor=tint(0.6), lw=0.7, zorder=2)
# highlight one real doublet and one T cell
for i, lab, dx, dy in ((500 + int(np.argmax(pann[500:])), None, 0, 0),
                       (int(np.argmin(pann[:300] + (np.linalg.norm(T - [-2.4, -0.3], axis=1) > 0.6))), None, 0, 0)):
    for j in idx[i, 1:]:
        b.plot([real[i, 0], allp[j, 0]], [real[i, 1], allp[j, 1]], color=TEXT2, lw=0.5, zorder=1)
    b.scatter([real[i, 0]], [real[i, 1]], s=110, facecolor="none", edgecolor=TEXT, lw=1, zorder=4)
i1 = 500 + int(np.argmax(pann[500:])); i2 = int(np.argmin(pann[:300] + (np.linalg.norm(T - [-2.4, -0.3], axis=1) > 0.6)))
b.annotate(f"pANN = {pann[i1]:.2f}\n→ doublet", real[i1], xytext=(-0.5, 2.5), fontsize=8.6, color=CRIMSON,
           fontweight="semibold", ha="center", arrowprops=dict(arrowstyle="-", color=CRIMSON, lw=0.7), linespacing=1.15)
b.annotate(f"pANN = {pann[i2]:.2f}\n→ singlet", real[i2], xytext=(-2.4, -2.6), fontsize=8.6, color=TEXT,
           ha="center", arrowprops=dict(arrowstyle="-", color=TEXT2, lw=0.7), linespacing=1.15)
b.text(-2.2, 1.35, "T cells", ha="center", fontsize=8.6, color=TEXT2); b.text(2.2, 2.0, "monocytes", ha="center", fontsize=8.6, color=TEXT2)
b.set_xlim(-4.2, 4.2); b.set_ylim(-3.4, 3.4); b.set_xticks([]); b.set_yticks([])
for sp in ("left", "bottom"): b.spines[sp].set_color(LINE)
b.set_xlabel("PC space", color=TEXT2, fontsize=8.8)
b.set_title("② Embed together, count neighbours", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
b.text(0, 1.03, f"pANN = share of artificial doublets among k = {k} neighbours", transform=b.transAxes, fontsize=9, color=TEXT2)
from matplotlib.lines import Line2D
b.legend([Line2D([], [], marker="o", ls="", color=CRIMSON, ms=4.5),
          Line2D([], [], marker="D", ls="", mfc="none", mec=tint(0.6), ms=4.5)],
         ["real doublets (hidden truth)", "artificial doublets"], loc="lower right", frameon=False, fontsize=8,
         handletextpad=0.2, labelcolor=TEXT2, borderaxespad=0.1)

# (c) rank by pANN and cut
c = fig.add_subplot(gs[2])
bins = np.linspace(0, 1, 17)
c.hist(pann[kind < 2], bins=bins, color=MC, lw=0, label="singlets")
c.hist(pann[kind == 2], bins=bins, color=CRIMSON, lw=0, label="doublets")
cut = np.sort(pann)[::-1][24]
c.axvline(cut + 0.03, color=TEXT, lw=0.9, ls=(0, (4, 3)))
c.text(cut + 0.06, c.get_ylim()[1] * 0.92, "cut: top n = expected\ndoublet rate × cells", fontsize=8.3,
       color=TEXT, va="top", linespacing=1.2)
c.set_yscale("log"); c.set_xlabel("pANN", color=TEXT2, fontsize=9); c.set_ylabel("cells", color=TEXT2, fontsize=9)
for sp in ("left", "bottom"): c.spines[sp].set_color(LINE)
c.tick_params(colors=TEXT2, labelsize=8); c.minorticks_off()
c.legend(frameon=False, fontsize=8.3, loc="center right", labelcolor=TEXT2, handlelength=1)
c.set_title("③ Rank by pANN and cut", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
c.text(0, 1.03, "two humps → separable", transform=c.transAxes, fontsize=9, color=TEXT2)
fig.text(0.02, 0.0, "A T cell + T cell doublet lands inside the T-cell cloud: homotypic doublets are invisible to this test.",
         fontsize=8.8, color=TEXT2)
header(a, "Single-cell RNA-seq · doublets", "Make fake doublets, then ask which real cells look like them",
       dy=40, x1=3.75)
fig.savefig(OUT + "doublet_detection.png"); print("ok k11 doublet", pann[kind == 2].mean().round(2), pann[kind < 2].mean().round(2))
