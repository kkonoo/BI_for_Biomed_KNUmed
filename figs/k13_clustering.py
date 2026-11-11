from gallery import *
from sklearn.neighbors import NearestNeighbors
rng = np.random.default_rng(5)

# Clustering = cutting the graph where links are sparse; Louvain vs Leiden; judge clusters by stability
W, H = 9.8, 4.5
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, (H - 0.8) / H]); ax.set_xlim(0, W); ax.set_ylim(0, H - 0.8); ax.axis("off")
DARK = "#5E6670"; MID = "#A9B0B9"; COL = [DARK, tint(0.5), MID]
def ptitle(x, t, sub):
    ax.text(x, 3.5, t, ha="left", va="center", fontsize=11.6, fontweight="semibold", color=TEXT)
    ax.text(x, 3.24, sub, ha="left", va="center", fontsize=8.6, color=TEXT2)

# (a) community detection on a small graph
ptitle(0.1, "(a) Community detection", "cut where links are sparse")
cen = np.array([[0.8, 2.45], [2.3, 2.5], [1.55, 1.35]])
X = np.vstack([rng.normal(c, 0.22, (16, 2)) for c in cen]); g = np.repeat([0, 1, 2], 16)
_, idx = NearestNeighbors(n_neighbors=5).fit(X).kneighbors(X)
E = {(min(i, j), max(i, j)) for i, r in enumerate(idx) for j in r[1:]}
E = {e for e in E if g[e[0]] == g[e[1]]}
for a_, b_ in ((0, 1), (1, 2), (0, 2)):                     # one or two sparse links between communities
    i = np.argmin(np.linalg.norm(X[g == a_] - cen[b_], axis=1)) + 16 * a_
    j = np.argmin(np.linalg.norm(X[g == b_] - cen[a_], axis=1)) + 16 * b_
    ax.plot(*zip(X[i], X[j]), color=CRIMSON, lw=1, ls=(0, (2, 1.5)), zorder=1)
for i, j in E: ax.plot(*zip(X[i], X[j]), color=COL[g[i]], lw=0.8, zorder=1)
ax.scatter(X[:, 0], X[:, 1], s=16, c=np.array(COL)[g], lw=0, zorder=2)
ax.text(0.1, 0.45, "more links inside than chance predicts\n(modularity); red = the links that get cut",
        ha="left", va="center", fontsize=8.3, color=TEXT2, linespacing=1.2)

# (b) Louvain vs Leiden
x0 = 3.35
ptitle(x0, "(b) Louvain vs Leiden", "Leiden adds a refinement step")
def mini(y, cols, verdict, good):
    P = {"l1": (0.15, 0.25), "l2": (0.0, -0.05), "l3": (0.3, -0.05), "M": (0.9, 0.35), "r1": (1.65, 0.25),
         "r2": (1.5, -0.05), "r3": (1.8, -0.05), "o1": (0.7, 0.75), "o2": (1.1, 0.75)}
    edges = [("l1", "l2"), ("l2", "l3"), ("l1", "l3"), ("r1", "r2"), ("r2", "r3"), ("r1", "r3"), ("l1", "M"),
             ("r1", "M"), ("M", "o1"), ("M", "o2"), ("o1", "o2")]
    for a_, b_ in edges:
        ax.plot([x0 + 0.3 + P[a_][0], x0 + 0.3 + P[b_][0]], [y + P[a_][1], y + P[b_][1]], color=LINE, lw=1, zorder=1)
    for n, (px, py) in P.items():
        ax.scatter([x0 + 0.3 + px], [y + py], s=60, color=cols[n], lw=0, zorder=2)
    ax.text(x0 + 2.35, y + 0.2, verdict, ha="left", va="center", fontsize=8.5, linespacing=1.2,
            color=TEXT if good else CRIMSON, fontweight="semibold")
louv = {**{k: DARK for k in ("l1", "l2", "l3", "r1", "r2", "r3")}, "M": tint(0.5), "o1": tint(0.5), "o2": tint(0.5)}
leid = {**{k: DARK for k in ("l1", "l2", "l3")}, **{k: MID for k in ("r1", "r2", "r3")}, "M": tint(0.5),
        "o1": tint(0.5), "o2": tint(0.5)}
ax.text(x0, 2.85, "Louvain", ha="left", va="center", fontsize=9.3, fontweight="bold", color=TEXT)
mini(1.95, louv, "✕ M moved away:\n  the grey cluster is\n  now two pieces", False)
ax.text(x0, 1.45, "Leiden", ha="left", va="center", fontsize=9.3, fontweight="bold", color=TEXT)
mini(0.55, leid, "✓ every cluster\n  is connected", True)
ax.text(x0, 0.12, "biological conclusions rarely change", ha="left", va="center", fontsize=8.3, color=TEXT2)

# (c) stability across resolutions (clustree-style)
x0 = 6.75
ptitle(x0, "(c) Check stability, not the number", "follow clusters as resolution rises")
rows = [("0.2", [("A", 600, 0.7), ("B", 400, 2.2)]),
        ("0.5", [("A1", 340, 0.35), ("A2", 260, 1.15), ("B", 400, 2.2)]),
        ("1.0", [("A1", 340, 0.35), ("A2", 200, 1.05), ("B1", 250, 1.85), ("B2", 210, 2.6)])]
pos = {}
for r, (res, nodes) in enumerate(rows):
    y = 2.75 - r * 0.95
    ax.text(x0, y, res, ha="left", va="center", fontsize=8.3, color=TEXT2)
    for nm, n, x in nodes:
        pos[(r, nm)] = (x0 + 0.4 + x * 0.95, y)
links = [((0, "A"), (1, "A1"), 340, 0), ((0, "A"), (1, "A2"), 260, 0), ((0, "B"), (1, "B"), 400, 0),
         ((1, "A1"), (2, "A1"), 340, 0), ((1, "A2"), (2, "A2"), 200, 0), ((1, "A2"), (2, "B1"), 60, 1),
         ((1, "B"), (2, "B1"), 190, 1), ((1, "B"), (2, "B2"), 210, 1)]
for a_, b_, n, bad in links:
    ax.plot(*zip(pos[a_], pos[b_]), color=CRIMSON if bad else LINE, lw=0.6 + n / 90, zorder=1, alpha=0.9)
for r, (res, nodes) in enumerate(rows):
    for nm, n, x in nodes:
        px, py = pos[(r, nm)]; bad = r == 2 and nm.startswith("B")
        ax.scatter([px], [py], s=n / 1.1, color=tint(0.55) if bad else DARK, lw=0, zorder=2)
        ax.text(px, py, nm, ha="center", va="center", fontsize=7.2, color="white", fontweight="bold", zorder=3)
ax.text(W - 0.1, 0.3, "A1/A2: clean split, kept at every step → check markers\nB1/B2: swap cells with A2 → unstable, likely noise",
        ha="right", va="center", fontsize=8.1, color=TEXT2, linespacing=1.25)
header(ax, "Single-cell RNA-seq · clustering", "Clustering = cutting the graph where links are sparse", dy=0)
fig.savefig(OUT + "clustering.png"); print("ok k13 clustering")
