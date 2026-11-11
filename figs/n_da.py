"""Differential-abundance stand-ins for the W7 figure. Each is a simplified version of one approach, not the package:
same data in, different assumptions, different calls. Truth is defined on absolute (cell) abundance."""
import numpy as np
from scipy import stats


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def simulate(rng, n=40, ntax=120, fold=5.0):
    base = rng.normal(0, 1.5, ntax)
    top = np.argsort(-base)
    up = np.r_[top[:3], rng.choice(top[10:], 5, replace=False)]
    rest = np.setdiff1d(np.arange(ntax), up)
    down = rng.choice(rest[rest > 0], 5, replace=False)
    eff = np.zeros(ntax); eff[up] = np.log(fold); eff[down] = -np.log(fold)
    g = np.r_[np.zeros(n, int), np.ones(n, int)]
    logc = base + rng.normal(0, 1.0, (2 * n, ntax)) + np.outer(g, eff)
    cells = np.exp(logc) * (rng.random((2 * n, ntax)) < 0.75)
    load = cells.sum(1)
    depth = np.exp(rng.normal(np.log(15000), 0.5, 2 * n)).astype(int)
    counts = np.array([rng.multinomial(d, c / c.sum()) for d, c in zip(depth, cells)])
    return dict(counts=counts, g=g, truth=eff != 0, load=load, up=up, down=down)


def methods(S, rng, alpha=0.05):
    C, g = S["counts"].astype(float), S["g"]
    P = C / C.sum(1, keepdims=True)
    a, b = g == 0, g == 1
    out = {}
    # 1 proportions + Wilcoxon, BH
    p1 = stats.mannwhitneyu(P[b], P[a], axis=0).pvalue
    out["proportions + Wilcoxon"] = bh(p1) < alpha
    # 2 proportions, Wilcoxon without FDR correction (as LEfSe's default)
    out["proportions, no FDR\n(LEfSe-style)"] = p1 < alpha
    # 3 median-of-ratios size factors + t-test on log counts (DESeq2-style normalisation)
    pos = np.where(C > 0, C, np.nan)
    gmean = np.exp(np.nanmean(np.log(pos), 0))
    sf = np.nanmedian(pos / gmean, 1)
    L = np.log2(C / sf[:, None] + 1)
    p3 = stats.ttest_ind(L[b], L[a], axis=0, equal_var=False).pvalue
    out["median-of-ratios + t\n(DESeq2-style)"] = bh(np.nan_to_num(p3, nan=1)) < alpha
    # 4 Monte Carlo Dirichlet + CLR + Wilcoxon, expected BH q (ALDEx2-style)
    qs = []
    for _ in range(32):
        D = np.array([rng.dirichlet(c + 0.5) for c in C])
        Z = np.log(D); Z -= Z.mean(1, keepdims=True)
        qs.append(bh(stats.mannwhitneyu(Z[b], Z[a], axis=0).pvalue))
    out["CLR + Wilcoxon\n(ALDEx2-style)"] = np.mean(qs, 0) < alpha
    # 5 CLR + linear model, shift corrected by the mode of all effects (LinDA-style)
    Z = np.log(C + 0.5); Z -= Z.mean(1, keepdims=True)
    d = Z[b].mean(0) - Z[a].mean(0)
    se = np.sqrt(Z[b].var(0, ddof=1) / b.sum() + Z[a].var(0, ddof=1) / a.sum())
    kde = stats.gaussian_kde(d); grid = np.linspace(d.min(), d.max(), 2001)
    mode = grid[np.argmax(kde(grid))]
    t = (d - mode) / se
    p5 = 2 * stats.t.sf(np.abs(t), len(g) - 2)
    out["CLR, shift corrected\n(LinDA-style)"] = bh(p5) < alpha
    # 6 reference: total load measured → absolute abundance
    A = np.log(P * S["load"][:, None] * np.exp(rng.normal(0, 0.15, len(g)))[:, None] + 1e-3 * S["load"].min())
    p6 = stats.ttest_ind(A[b], A[a], axis=0, equal_var=False).pvalue
    out["total load measured"] = bh(np.nan_to_num(p6, nan=1)) < alpha
    return out
