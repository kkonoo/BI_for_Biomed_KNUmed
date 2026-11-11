from gallery import *
from k_toy import *
import igraph as ig, random
random.seed(1); ig.set_random_number_generator(random)

# The cell-type workflow on the toy data: neighbour graph → clusters on the graph → UMAP picture → names
X, lab, lin, t = toy(scale=0.45); P, _ = pca(X)
idx = knn(P, 10); E = snn_edges(idx)
g = ig.Graph(n=len(P), edges=list(E)); g.es["weight"] = [w + 0.02 for w in E.values()]
lay = np.array(g.layout_fruchterman_reingold(weights="weight", niter=1500))
cl = leiden(P, 0.05, k=10); U = umap2(P)
fig, axes = plt.subplots(1, 4, figsize=(9.8, 3.4), gridspec_kw={"wspace": 0.06})
import collections
ccol = {c: LINCOL[collections.Counter(lin[cl == c]).most_common(1)[0][0]] for c in range(cl.max() + 1)}
cc = np.array([ccol[c] for c in cl])
def deco(ax, t, sub, hl=False):
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(True); sp.set_color("#E3E6EB")
    ax.set_title(t, loc="left", fontsize=11.3, fontweight="semibold", color=CRIMSON if hl else TEXT, pad=22)
    ax.text(0, 1.03, sub, transform=ax.transAxes, fontsize=8.4, color=TEXT2)
a, b, c, d = axes
for (i, j) in E:
    for ax in (a, b): ax.plot(*zip(lay[i], lay[j]), color="#D5D9DF", lw=0.3, zorder=1)
a.scatter(lay[:, 0], lay[:, 1], s=5, color="#5E6670", lw=0, zorder=2)
deco(a, "① Neighbour graph", "who is near whom in PC space")
b.scatter(lay[:, 0], lay[:, 1], s=5, c=cc, lw=0, zorder=2)
deco(b, "② Clustering", "cut the graph into communities", hl=True)
c.scatter(U[:, 0], U[:, 1], s=5, c=cc, lw=0)
deco(c, "③ UMAP", "a 2-D picture of the same graph")
d.scatter(U[:, 0], U[:, 1], s=5, c=cc, lw=0, alpha=0.5)
names = {0: "T / NK", 1: "B / plasma", 2: "DC", 3: "mono → mac", 4: "stromal"}
for c_ in range(cl.max() + 1):
    l = collections.Counter(lin[cl == c_]).most_common(1)[0][0]
    cx, cy = np.median(U[cl == c_], 0)
    d.text(cx, cy, names[l], ha="center", va="center", fontsize=8.3, fontweight="semibold", color=TEXT,
           bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.9))
deco(d, "④ Annotation", "names from marker genes")
for x in (0.305, 0.51, 0.715):
    fig.text(x, 0.45, "→", ha="center", va="center", fontsize=13, color=LINE)
fig.text(0.125, -0.04, "The graph has no real coordinates: ① and ② are one drawing of it.  Clustering runs on the graph, "
         "never on the UMAP picture.", fontsize=8.7, color=TEXT2)
header(a, "Single-cell RNA-seq · finding cell types", "Graph first, clusters on the graph, UMAP only to look",
       dy=40, x1=4.3)
fig.savefig(OUT + "celltype_workflow.png"); print("ok k24 flow", cl.max() + 1)
