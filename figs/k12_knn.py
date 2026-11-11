from gallery import *
from sklearn.neighbors import NearestNeighbors
rng = np.random.default_rng(12)

# Building the neighbour graph: k nearest neighbours in PC space, then weight edges by shared neighbours (SNN)
fig, axes = plt.subplots(1, 3, figsize=(9.8, 4.0), gridspec_kw={"wspace": 0.08})
DARK = "#5E6670"; GREYF = "#9AA1AB"
A = rng.normal([-1.5, 0.9], [0.33, 0.3], (22, 2))
t = rng.uniform(0, 1, 26); B = np.c_[-0.6 + 1.7 * t, 0.2 - 0.9 * t] + rng.normal(0, 0.16, (26, 2))
C = rng.normal([1.75, -1.15], [0.3, 0.28], (18, 2))
X = np.vstack([A, B, C]); k = 5
_, idx = NearestNeighbors(n_neighbors=k + 1).fit(X).kneighbors(X); idx = idx[:, 1:]
N = [set(r) | {i} for i, r in enumerate(idx)]; E = {}          # closed neighbourhoods; SNN weight = Jaccard
for i, r in enumerate(idx):
    for j in r:
        a, b = min(i, j), max(i, j); E[(a, b)] = len(N[a] & N[b]) / len(N[a] | N[b])
def base(ax, title, sub):
    ax.set_xlim(-2.5, 2.6); ax.set_ylim(-2.3, 1.85); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_color("#E3E6EB")
    ax.set_title(title, loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
    ax.text(0, 1.03, sub, transform=ax.transAxes, fontsize=8.9, color=TEXT2)
# ① one cell and its k neighbours
a = axes[0]; base(a, f"① Find each cell's k = {k}", "nearest cells in PC space")
a.scatter(X[:, 0], X[:, 1], s=16, color=GREYF, lw=0, zorder=2)
q = int(np.argmin(np.linalg.norm(X - [-1.5, 0.9], axis=1)))
r = np.linalg.norm(X[idx[q, -1]] - X[q])
a.add_patch(mp.Circle(X[q], r * 1.08, fc=tint(0.08), ec=CRIMSON, lw=0.9, ls=(0, (2, 2)), zorder=1))
for j in idx[q]: a.plot(*zip(X[q], X[j]), color=CRIMSON, lw=1, zorder=1)
a.scatter(X[idx[q], 0], X[idx[q], 1], s=22, color=tint(0.6), lw=0, zorder=3)
a.scatter(*X[q], s=34, color=CRIMSON, lw=0, zorder=4)
a.text(0.03, 0.03, "distance = Euclidean on the top PCs\n(Seurat default k = 20)", transform=a.transAxes, fontsize=8.2,
       color=TEXT2, linespacing=1.2)
# ② kNN graph
b = axes[1]; base(b, "② Link every cell to its k", "cells = nodes, neighbours = edges")
for (i, j) in E: b.plot(*zip(X[i], X[j]), color=GREYF, lw=0.8, zorder=1)
b.scatter(X[:, 0], X[:, 1], s=16, color=DARK, lw=0, zorder=2)
b.text(0.03, 0.03, "dense regions = candidate populations;\nshape is kept (elongated, branched)", transform=b.transAxes,
       fontsize=8.2, color=TEXT2, linespacing=1.2)
# ③ SNN weights
c = axes[2]; base(c, "③ Weight by shared neighbours", "SNN: strong if they share most neighbours")
for (i, j), w in E.items():
    weak = w <= 0.2
    c.plot(*zip(X[i], X[j]), color=CRIMSON if weak else DARK, lw=0.6 if weak else 0.4 + 2.6 * w,
           ls=(0, (2, 1.5)) if weak else "-", zorder=1)
c.scatter(X[:, 0], X[:, 1], s=16, color=DARK, lw=0, zorder=2)
c.text(0.03, 0.03, "links between dense regions share few\nneighbours → weak (red, dashed)", transform=c.transAxes,
       fontsize=8.2, color=CRIMSON, linespacing=1.2)
fig.text(0.5, -0.04, "k sets the scale: small k → fine structure, more fragments; large k → smoother graph, small groups absorbed",
         ha="center", fontsize=8.9, color=TEXT2)
header(axes[0], "Single-cell RNA-seq · neighbour graph", "Clustering and UMAP both start from this graph, not from the genes",
       dy=40, x1=3.17)
fig.savefig(OUT + "knn_graph.png"); print("ok k12 knn")
