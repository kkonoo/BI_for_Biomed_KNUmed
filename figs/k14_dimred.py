from gallery import *
from k_toy import *
from sklearn.manifold import TSNE

# PCA, t-SNE and UMAP of the same 50-dimensional toy cells
X, lab, lin, t = toy(); P, ev = pca(X)
T = TSNE(perplexity=30, init="random", random_state=1).fit_transform(P)
U = umap2(P)
fig, axes = plt.subplots(1, 3, figsize=(9.8, 4.9), gridspec_kw={"wspace": 0.1})
sub = {"CD8 T": 0.4, "NK": -0.45, "plasma": -0.4, "endothelial": -0.35}
cols = np.array([shade(LINCOL[l], sub.get(n, 0.001)) if n in sub else LINCOL[l] for n, l in zip(lab, lin)])
panels = [("PCA", "a linear shadow (PC1, PC2)", P[:, :2],
           ["✓ axes and distances are honest", "✕ sub-types overlap: 2 PCs cannot", "   hold 50 dimensions",
            "→ the input to every later step"]),
          ("t-SNE", "matches neighbour probabilities", T,
           ["✓ local neighbours kept", "✕ where each island lands is", "   arbitrary; slow", ""]),
          ("UMAP", "lays out the kNN graph", U,
           ["✓ tight islands, clear gaps;", "   the path stays a strand; fast", "✕ still a picture, not a", "   measurement"])]
for ax, (nm, s, E, notes) in zip(axes, panels):
    o = np.random.default_rng(0).permutation(len(E))
    ax.scatter(E[o, 0], E[o, 1], s=4, c=cols[o], lw=0, alpha=0.9, rasterized=True)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_color("#E3E6EB")
    ax.set_title(nm, loc="left", fontsize=12.5, fontweight="semibold", color=TEXT, pad=22)
    ax.text(0, 1.03, s, transform=ax.transAxes, fontsize=9, color=TEXT2)
    for i, n in enumerate(notes):
        ax.text(0, -0.06 - i * 0.058, n, transform=ax.transAxes, fontsize=8.6, va="top",
                color=CRIMSON if n.startswith("✕") or (i and notes[i - 1].startswith("✕") and n.startswith("   ")) else TEXT)
from matplotlib.lines import Line2D
hs = [Line2D([], [], marker="o", ls="", color=c, ms=5) for c in LINCOL]
fig.legend(hs, LINEAGES, loc="lower center", ncol=5, frameon=False, fontsize=8.5, bbox_to_anchor=(0.5, -0.24),
           labelcolor=TEXT, handletextpad=0.2, columnspacing=1.2)
header(axes[0], "Single-cell RNA-seq · embeddings", "Three ways to flatten the same 50-dimensional cells",
       dy=40, x1=3.2)
fig.savefig(OUT + "dimred_compare.png"); print("ok k14 dimred")
