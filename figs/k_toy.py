# Shared toy single-cell dataset for the W5 figures: 50 "genes", five lineages with sub-types,
# and one continuous monocyte → macrophage path. Plus PCA, kNN/SNN graph and Leiden helpers.
import numpy as np

LINEAGES = ["T / NK", "B / plasma", "DC", "monocyte → macrophage", "stromal"]

def toy(seed=0, scale=1.0):
    rng = np.random.default_rng(seed); L = 10
    major = {k: rng.normal(0, 1.6, L) for k in LINEAGES}
    spec = [("CD4 T", 0, 260, 0.35), ("CD8 T", 0, 200, 0.35), ("NK", 0, 120, 0.6), ("B", 1, 170, 0.4),
            ("plasma", 1, 60, 0.8), ("DC", 2, 45, 0.0), ("fibroblast", 4, 150, 0.5), ("endothelial", 4, 110, 0.5)]
    Z, lab, lin, t = [], [], [], []
    for nm, li, n, d in spec:
        n = int(n * scale); c = major[LINEAGES[li]] + rng.normal(0, d, L)
        Z.append(c + rng.normal(0, 0.5, (n, L))); lab += [nm] * n; lin += [li] * n; t += [np.nan] * n
    n = int(280 * scale); tt = rng.uniform(0, 1, n)
    a = major[LINEAGES[3]]; b = a + rng.normal(0, 2.6, L)
    Z.append(a + np.outer(tt, b - a) + rng.normal(0, 0.45, (n, L)))
    lab += ["monocyte"] * n; lin += [3] * n; t += list(tt)
    Z = np.vstack(Z); W = rng.normal(0, 1, (L, 50)) / np.sqrt(L)
    X = Z @ W + rng.normal(0, 0.9, (len(Z), 50))
    return X, np.array(lab), np.array(lin), np.array(t)

def pca(X, k=20):
    Xc = X - X.mean(0); U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    return U[:, :k] * S[:k], S ** 2 / (S ** 2).sum()

def knn(P, k=15):
    from sklearn.neighbors import NearestNeighbors
    _, idx = NearestNeighbors(n_neighbors=k + 1).fit(P).kneighbors(P)
    return idx[:, 1:]

def snn_edges(idx):
    """Edges of the kNN graph weighted by the Jaccard overlap of the two neighbour sets (SNN)."""
    sets = [set(r) for r in idx]; E = {}
    for i, r in enumerate(idx):
        for j in r:
            a, b = min(i, j), max(i, j)
            if (a, b) not in E:
                E[(a, b)] = len(sets[a] & sets[b]) / len(sets[a] | sets[b])
    return E

def leiden(P, res=1.0, k=15, seed=0):
    import igraph as ig, leidenalg as la
    E = snn_edges(knn(P, k)); g = ig.Graph(n=len(P), edges=list(E)); g.es["weight"] = list(E.values())
    part = la.find_partition(g, la.RBConfigurationVertexPartition, weights="weight", resolution_parameter=res, seed=seed)
    return np.array(part.membership)

def umap2(P, seed=0, **kw):
    import umap, warnings
    warnings.filterwarnings("ignore")
    kw = {"n_neighbors": 15, "min_dist": 0.3, **kw}
    return umap.UMAP(random_state=seed, **kw).fit_transform(P)

# lineage colours in the gallery palette (crimson family + greys) and a shade helper for sub-groups
LINCOL = ["#5E6670", "#D88D99", "#2B2F36", "#B0283C", "#A9B0B9"]

def shade(c, f):
    """f > 0 mixes towards white, f < 0 towards black."""
    c = np.array([int(c[i:i + 2], 16) for i in (1, 3, 5)], float)
    t = np.full(3, 255.0) if f > 0 else np.zeros(3)
    return "#%02X%02X%02X" % tuple((c + (t - c) * abs(f)).round().astype(int))
