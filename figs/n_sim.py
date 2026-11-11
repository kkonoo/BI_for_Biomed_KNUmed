"""Shared helpers for the W7 figures: a simple gut-community simulator, Bray–Curtis, PCoA and PERMANOVA."""
import numpy as np
from scipy.special import gammaln


def _draw(n, rng, base, effect, n_eff, person_sd, carry, depth):
    ntaxa = len(base)
    logp = base + rng.normal(0, person_sd, (n, ntaxa))
    logp[:, :n_eff] += effect
    present = rng.random((n, ntaxa)) < carry
    p = np.exp(logp) * present
    p /= p.sum(1, keepdims=True)
    return np.array([rng.multinomial(depth, q) for q in p])


def two_groups(n_per, rng, effect, **kw):
    """Controls and cases drawn from the same base community; cases have the first n_eff taxa shifted."""
    ntaxa = kw.pop("ntaxa", 150)
    base = rng.normal(0, kw.get("base_sd", 1.2), ntaxa)
    args = dict(n_eff=kw.get("n_eff", 10), person_sd=kw.get("person_sd", 1.2), carry=kw.get("carry", 0.7),
                depth=kw.get("depth", 20000))
    c0 = _draw(n_per, rng, base, 0.0, **args)
    c1 = _draw(n_per, rng, base, effect, **args)
    return np.vstack([c0, c1]), np.r_[np.zeros(n_per, int), np.ones(n_per, int)]


def bray_curtis(C):
    P = C / C.sum(1, keepdims=True)
    n = len(P); D = np.zeros((n, n))
    for i in range(n):
        D[i] = 0.5 * np.abs(P - P[i]).sum(1)
    return D


def pcoa(D, k=2):
    n = len(D); J = np.eye(n) - 1 / n
    B = -0.5 * J @ (D ** 2) @ J
    w, v = np.linalg.eigh(B)
    o = np.argsort(w)[::-1]; w, v = w[o], v[:, o]
    expl = w[:k] / w[w > 0].sum()
    return v[:, :k] * np.sqrt(w[:k]), expl


def permanova(D, g, rng, nperm=999):
    """Two-group PERMANOVA (Anderson 2001). Returns pseudo-F, R², permutation p."""
    n = len(g); D2 = D ** 2
    sst = D2.sum() / (2 * n)

    def ssw(L):  # L: (m × n) 0/1 labels
        out = 0
        for k in (0, 1):
            V = (L == k).astype(float)
            out = out + np.einsum("ij,ij->i", V @ D2, V) / (2 * V.sum(1))
        return out
    w_obs = ssw(g[None, :])[0]
    F = lambda w: ((sst - w) / 1) / (w / (n - 2))
    perms = np.array([rng.permutation(g) for _ in range(nperm)])
    Fp = F(ssw(perms))
    f_obs = F(w_obs)
    return f_obs, (sst - w_obs) / sst, (np.sum(Fp >= f_obs) + 1) / (nperm + 1)


def rarefied_richness(counts, depths):
    """Expected number of taxa seen when subsampling `d` reads without replacement (Hurlbert 1971)."""
    N = counts.sum(); c = counts[counts > 0]
    lc = lambda a, b: gammaln(a + 1) - gammaln(b + 1) - gammaln(a - b + 1)
    out = []
    for d in depths:
        ok = (N - c) >= d
        miss = np.zeros_like(c, float)
        miss[ok] = np.exp(lc(N - c[ok], d) - lc(N, d))
        out.append((1 - miss).sum())
    return np.array(out)
