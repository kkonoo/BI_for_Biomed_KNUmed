from gallery import *
rng = np.random.default_rng(23)

# Normalise → highly variable genes → scale → PCA, computed on simulated UMI counts (5 cell types, 3,000 genes)
nc, ng, nt = 400, 3000, 5
typ = rng.integers(0, nt, nc)
base = np.exp(rng.normal(-1.3, 1.5, ng))
fold = np.ones((nt, ng)); mk = rng.choice(ng, 250, replace=False)
for i, g in enumerate(mk): fold[i % nt, g] = rng.uniform(4, 12)
sf = np.exp(rng.normal(0, 0.5, nc))
mu = sf[:, None] * base[None] * fold[typ]
counts = rng.poisson(rng.gamma(1 / 0.15, mu * 0.15))                       # negative binomial, dispersion 0.15
tot = counts.sum(1)
L = np.log1p(counts / tot[:, None] * 1e4)
m, v = L.mean(0), L.var(0)
ok = m > 0.02
o = np.argsort(m[ok]); bins = np.array_split(o, 20)                         # variance relative to genes of similar mean
z = np.zeros(ok.sum())
for bb in bins: z[bb] = (v[ok][bb] - v[ok][bb].mean()) / (v[ok][bb].std() + 1e-9)
hv = np.zeros(ng, bool); hv[np.where(ok)[0][np.argsort(-z)[:300]]] = True
S = L[:, hv]; S = (S - S.mean(0)) / S.std(0)
U_, s_, Vt = np.linalg.svd(S, full_matrices=False); ev = s_ ** 2 / (s_ ** 2).sum()
DARK = "#5E6670"; GREYF = "#C3C8CF"
fig, axes = plt.subplots(1, 4, figsize=(9.8, 3.7), gridspec_kw={"wspace": 0.45, "width_ratios": [1, 1, 1, 1.05]})
def deco(ax, t, sub):
    for sp in ("left", "bottom"): ax.spines[sp].set_color(LINE)
    ax.tick_params(colors=TEXT2, labelsize=7.8)
    ax.set_title(t, loc="left", fontsize=11.2, fontweight="semibold", color=TEXT, pad=22)
    ax.text(0, 1.03, sub, transform=ax.transAxes, fontsize=8.3, color=TEXT2)
# ① normalise
a = axes[0]; k = rng.choice(nc, 8, replace=False)
a.bar(np.arange(8) - 0.2, tot[k], width=0.4, color=GREYF, label="raw")
a.bar(np.arange(8) + 0.2, np.full(8, 1e4), width=0.4, color=CRIMSON, label="normalised")
a.set_xticks([]); a.set_xlabel("8 cells", color=TEXT2, fontsize=8.5); a.set_ylabel("total UMI", color=TEXT2, fontsize=8.5)
a.legend(frameon=False, fontsize=7.6, loc="upper right", labelcolor=TEXT2, handlelength=0.8, bbox_to_anchor=(1.05, 1.02))
a.set_ylim(0, 1.45e4)
deco(a, "① Normalise", "÷ cell total × 10⁴, log(1+x)")
a.text(0, -0.2, f"raw totals differ {tot[k].max() / tot[k].min():.1f}×,\nfor technical reasons", transform=a.transAxes,
       fontsize=8, color=TEXT2, va="top", linespacing=1.2)
# ② HVG
b = axes[1]
b.scatter(m[ok & ~hv], v[ok & ~hv], s=3, color=GREYF, lw=0, rasterized=True)
b.scatter(m[hv], v[hv], s=5, color=CRIMSON, lw=0, rasterized=True)
b.set_xscale("log"); b.set_yscale("log"); b.minorticks_off()
b.set_xlabel("mean", color=TEXT2, fontsize=8.5); b.set_ylabel("variance", color=TEXT2, fontsize=8.5)
deco(b, "② Keep variable genes", "variance high for its mean")
b.text(0, -0.2, f"red = {hv.sum()} of {ng:,} genes here\n(~2,000 of 20,000 in practice)", transform=b.transAxes,
       fontsize=8, color=TEXT2, va="top", linespacing=1.2)
# ③ scale + PCA: elbow
c = axes[2]; n = 20
c.plot(np.arange(1, n + 1), 100 * ev[:n], "o-", color=DARK, ms=3.5, lw=1.2)
c.plot(np.arange(1, nt), 100 * ev[:nt - 1], "o", color=CRIMSON, ms=4.5)
c.set_xlabel("PC", color=TEXT2, fontsize=8.5); c.set_ylabel("% variance", color=TEXT2, fontsize=8.5)
c.set_xticks([1, 5, 10, 15, 20])
deco(c, "③ Scale, then PCA", "each gene mean 0, SD 1")
c.text(0, -0.2, f"first {nt - 1} PCs (red) = the {nt} cell types;\nthe flat tail is mostly noise", transform=c.transAxes,
       fontsize=8, color=TEXT2, va="top", linespacing=1.2)
# ④ a PC is a gene program
d = axes[3]
ld = Vt[0] * np.sign(Vt[0][np.argmax(np.abs(Vt[0]))]); idx = np.argsort(-np.abs(ld))[:10][::-1]
gid = np.where(hv)[0][idx]
ismk = np.isin(gid, mk)
d.barh(np.arange(10), ld[idx], color=[CRIMSON if x else GREYF for x in ismk], height=0.7)
d.set_yticks([]); d.set_ylabel("top 10 genes on PC1", color=TEXT2, fontsize=8.5)
d.axvline(0, color=LINE, lw=0.8); d.set_xlabel("PC1 loading", color=TEXT2, fontsize=8.5)
deco(d, "④ A PC = a gene program", "genes that move together")
d.text(0, -0.2, f"{ismk.sum()} of 10 are marker genes\nof one cell type (red)", transform=d.transAxes, fontsize=8,
       color=TEXT2, va="top", linespacing=1.2)
fig.text(0.06, -0.21, "Graph, clustering and UMAP all run on 10–50 PCs.  Tests go back to the raw counts (§6).",
         fontsize=8.8, color=TEXT)
header(axes[0], "Single-cell RNA-seq · normalisation and PCA", "From 20,000 noisy genes to a few dozen stable axes",
       dy=40, x1=5.6)
fig.savefig(OUT + "norm_pca.png"); print("ok k23 norm pca", ev[:6].round(3), ismk.sum())
