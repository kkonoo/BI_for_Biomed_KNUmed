# Shared simulation for the W6 missing-value figures (m6, m12).
# 3,000 proteins, 6 vs 6 samples, log2 intensities. 20% of proteins change, by at most 0.8 log2 (1.7-fold).
# A value goes missing if it falls below the detection limit (MNAR, soft cut-off at 21)
# or, at any intensity, by chance (MAR, 5%).
import numpy as np

def missing_sim(seed=0, P=3000, n=6, lod=21.0):
    rng = np.random.default_rng(seed)
    mu = rng.normal(23, 2.2, P)
    tl = np.where(rng.random(P) < 0.2, rng.choice([-1, 1], P) * rng.uniform(0.2, 0.8, P), 0)
    X = mu[:, None] + np.c_[np.zeros((P, n)), np.repeat(tl[:, None], n, 1)] + rng.normal(0, 0.4, (P, 2 * n))
    below = rng.random(X.shape) > 1 / (1 + np.exp(-(X - lod) / 0.4))
    chance = rng.random(X.shape) < 0.05
    O = np.where(below | chance, np.nan, X)
    return dict(X=X, O=O, tl=tl, n=n, below=below, chance=chance & ~below, rng=rng, lod=lod)

def bh(p):
    p = np.asarray(p, float); m = len(p); o = np.argsort(p)
    q = p[o] * m / np.arange(1, m + 1); q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(m); out[o] = np.minimum(q, 1); return out
