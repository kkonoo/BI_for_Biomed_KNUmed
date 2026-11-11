from gallery import *
from k_toy import *

# Reading a UMAP: neighbourhoods are trustworthy; distances between islands, island sizes and bridges are not.
# Toy data + 30 artificial T/NK–B doublets; distances and spreads are measured, not drawn.
X, lab, lin, t = toy(); rng = np.random.default_rng(25)
T_, B_ = np.where(lin == 0)[0], np.where(lin == 1)[0]
w = rng.uniform(0.1, 0.9, (40, 1))                       # doublets with uneven mixing span the gap
D = w * X[rng.choice(T_, 40)] + (1 - w) * X[rng.choice(B_, 40)]
X2 = np.vstack([X, D]); lin2 = np.r_[lin, np.full(40, 9)]
P, _ = pca(X2); U = umap2(P)
C = {l: P[lin2 == l].mean(0) for l in range(5)}; CU = {l: np.median(U[lin2 == l], 0) for l in range(5)}
dpc = lambda a, b: np.linalg.norm(C[a] - C[b]); dum = lambda a, b: np.linalg.norm(CU[a] - CU[b])
spc = lambda l: np.linalg.norm(P[lin2 == l] - C[l], axis=1).mean()
sum_ = lambda l: np.linalg.norm(U[lin2 == l] - CU[l], axis=1).mean()
nnP = knn(P, 15); nnU = knn(U, 15)
keep = np.mean([len(set(nnP[i]) & set(nnU[i])) / 15 for i in range(len(X))])

fig = plt.figure(figsize=(9.8, 5.0))
gs = fig.add_gridspec(2, 2, width_ratios=[1.25, 1], height_ratios=[1, 1], wspace=0.28, hspace=0.95,
                      left=0.01, right=0.98, top=0.8, bottom=0.1)
a = fig.add_subplot(gs[:, 0])
cols = np.array([LINCOL[l] if l < 9 else CRIMSON for l in lin2])
a.scatter(U[:, 0], U[:, 1], s=4, c=cols, lw=0, alpha=0.55, rasterized=True)
a.scatter(U[lin2 == 9, 0], U[lin2 == 9, 1], s=12, color=CRIMSON, lw=0, zorder=3)
q = int(np.where(lin2 == 0)[0][5])
a.scatter(U[nnP[q], 0], U[nnP[q], 1], s=26, facecolor="none", edgecolor=TEXT, lw=0.8, zorder=4)
a.scatter(*U[q], s=40, color=TEXT, zorder=5)
for l, nm in enumerate(["T / NK", "B / plasma", "DC", "mono → mac", "stromal"]):
    if l == 3:
        a.text(U[lin2 == l, 0].min() - 0.3, CU[l][1], nm, ha="right", va="center", fontsize=8.3, color=TEXT2)
    else:
        a.text(CU[l][0], U[lin2 == l, 1].max() + 0.35, nm, ha="center", va="bottom", fontsize=8.3, color=TEXT2)
bm = np.median(U[lin2 == 9], 0)
a.annotate("✕ a bridge made of\n   40 doublets", bm, xytext=(0.57, 0.6), textcoords="axes fraction", linespacing=1.15, fontsize=8.5,
           color=CRIMSON, fontweight="semibold", arrowprops=dict(arrowstyle="-", color=CRIMSON, lw=0.7))
a.annotate("✓ one cell (●) and its 15 nearest\n   neighbours in PC space (○)\n   land close together",
           U[q], xytext=(0.55, 0.22), textcoords="axes fraction", fontsize=8.5, color=TEXT, linespacing=1.15,
           arrowprops=dict(arrowstyle="-", color=TEXT2, lw=0.7))
a.set_xticks([]); a.set_yticks([])
for sp in a.spines.values(): sp.set_visible(True); sp.set_color("#E3E6EB")
a.set_title("What you can trust: who sits next to whom", loc="left", fontsize=12, fontweight="semibold", color=TEXT, pad=22)
a.text(0, 1.03, "toy data (50 genes) with 40 T/NK + B doublets added", transform=a.transAxes, fontsize=8.8, color=TEXT2)

def pair(ax, title, sub, labels, v_pc, v_um, note):
    v_pc = np.array(v_pc) / max(v_pc); v_um = np.array(v_um) / max(v_um)
    y = np.arange(len(labels))[::-1]
    ax.barh(y + 0.18, v_pc, height=0.34, color="#5E6670", label="in PC space")
    ax.barh(y - 0.18, v_um, height=0.34, color=CRIMSON, label="on the UMAP")
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8.4, color=TEXT)
    ax.set_xticks([]); ax.set_xlim(0, 1.25)
    for sp in ("left",): ax.spines[sp].set_color(LINE)
    ax.spines["bottom"].set_visible(False); ax.tick_params(axis="y", length=0)
    ax.set_title(title, loc="left", fontsize=11.5, fontweight="semibold", color=CRIMSON, pad=20)
    ax.text(0, 1.04, sub, transform=ax.transAxes, fontsize=8.5, color=TEXT2)
    ax.text(0, -0.2, note, transform=ax.transAxes, fontsize=8.2, color=TEXT2, va="top")
b = fig.add_subplot(gs[0, 1])
pair(b, "✕ Distances between islands", "grey = in PC space, red = on the UMAP",
     ["T / NK ↔ stromal", "B / plasma ↔ stromal"], [dpc(0, 4), dpc(1, 4)], [dum(0, 4), dum(1, 4)],
     "distance between centres, each scaled to the larger:\nthe order flips between PC space and the UMAP")
c = fig.add_subplot(gs[1, 1])
pair(c, "✕ Island sizes", "mean distance of cells to their island centre",
     ["T / NK", "DC"], [spc(0), spc(2)], [sum_(0), sum_(2)],
     f"DC: {spc(2) / spc(0):.1f}× the spread of T / NK in PC space, {sum_(2) / sum_(0):.2f}× on the UMAP")
header(a, "Single-cell RNA-seq · reading a UMAP", "Trust the neighbourhoods; do not measure the map", dy=40, x1=1.85)
fig.savefig(OUT + "umap_reading.png")
print("ok k25", round(keep, 2), round(dpc(0, 4), 1), round(dpc(1, 4), 1), round(dum(0, 4), 1), round(dum(1, 4), 1),
      round(spc(2) / spc(0), 2), round(sum_(2) / sum_(0), 2))
