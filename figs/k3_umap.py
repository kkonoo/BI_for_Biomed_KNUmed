from gallery import *
from k_toy import *
import collections

# Same cells, same algorithm, one parameter moved: Leiden at a low and a high resolution, shown on one UMAP
X, lab, lin, t = toy(); P, _ = pca(X); U = umap2(P)
lo = leiden(P, 0.05); hi = leiden(P, 2.0)
fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.9), gridspec_kw={"wspace": 0.06})
# low-res clusters take their lineage colour; each high-res cluster takes a shade of its parent's colour
pcol = {c: LINCOL[collections.Counter(lin[lo == c]).most_common(1)[0][0]] for c in range(lo.max() + 1)}
kids = collections.defaultdict(list)
for c in range(hi.max() + 1):
    kids[collections.Counter(lo[hi == c]).most_common(1)[0][0]].append(c)
hcol = {}
for p, cs in kids.items():
    cs = sorted(cs, key=lambda c: -np.sum(hi == c))
    for i, c in enumerate(cs):
        hcol[c] = pcol[p] if i == 0 else shade(pcol[p], [-0.4, 0.5, -0.65, 0.75][(i - 1) % 4])
for ax, m, col, nm, sub in ((axes[0], lo, pcol, "Low resolution", "one cluster per lineage"),
                            (axes[1], hi, hcol, "High resolution", "lineages split; the continuous path is chopped into pieces")):
    for c in range(m.max() + 1):
        s = m == c
        ax.scatter(U[s, 0], U[s, 1], s=5, color=col[c], lw=0, alpha=0.9, rasterized=True)
        cx, cy = np.median(U[s], 0)
        ax.text(cx, cy, str(c + 1), ha="center", va="center", fontsize=8.5, fontweight="bold", color=TEXT,
                bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.85))
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.set_title(f"{nm} → {m.max() + 1} clusters", loc="left", fontsize=12.5, fontweight="semibold",
                 color=CRIMSON if m is hi else TEXT, pad=22)
    ax.text(0, 1.03, sub, transform=ax.transAxes, fontsize=9, color=TEXT2)
fig.text(0.5, 0.0, "No statistic says which is right: check marker genes and stability, and report the value you used",
         ha="center", fontsize=9, color=TEXT2)
header(axes[0], "Single-cell RNA-seq · clustering", "Same cells, same algorithm, one parameter moved", dy=40, x1=2.06)
fig.savefig(OUT + "umap_resolution.png"); print("ok k3 umap", lo.max() + 1, hi.max() + 1)
