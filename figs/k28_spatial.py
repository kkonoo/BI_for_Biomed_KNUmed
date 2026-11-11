from gallery import *
rng = np.random.default_rng(28)

# Spatial transcriptomics: spot-based (whole transcriptome, several cells per spot) vs imaging-based
# (single-cell or subcellular, a fixed gene panel)
fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.6), gridspec_kw={"wspace": 0.1})
DARK = "#5E6670"; GREYF = "#B8BEC7"; TC = [DARK, tint(0.5), GREYF]
gx, gy = np.meshgrid(np.arange(0, 8.4, 0.42), np.arange(0, 5.9, 0.42))
cells = np.c_[gx.ravel() + rng.normal(0, 0.06, gx.size), gy.ravel() + rng.normal(0, 0.06, gy.size)]
cells[:, 0] += (np.round(cells[:, 1] / 0.42) % 2) * 0.21
ct = np.where(cells[:, 1] > 3.9 + 0.4 * np.sin(cells[:, 0]), 0, np.where(cells[:, 1] > 1.9, 1, 2))
for ax in axes:
    for (x, y), t in zip(cells, ct):
        ax.add_patch(mp.Circle((x, y), 0.18, fc=TC[t], ec="white", lw=0.6, alpha=0.9))
    ax.set_xlim(-0.4, 8.7); ax.set_ylim(-1.4, 6.0); ax.set_aspect("equal"); ax.axis("off")
# (a) spots
a = axes[0]
spots = [(x, y) for i, y in enumerate(np.arange(0.5, 5.8, 1.25)) for x in np.arange(0.5 + (i % 2) * 0.72, 8.0, 1.44)]
for x, y in spots: a.add_patch(mp.Circle((x, y), 0.55, fc="none", ec=CRIMSON, lw=1.1))
def ent(p):                                              # pick the spot whose cells are the most mixed
    m = np.bincount(ct[np.linalg.norm(cells - p, axis=1) < 0.6], minlength=3); m = m[m > 0] / m.sum()
    return -(m * np.log(m)).sum()
hx, hy = max(spots, key=ent)
a.add_patch(mp.Circle((hx, hy), 0.55, fc=tint(0.12), ec=CRIMSON, lw=2, alpha=0.6))
inside = np.linalg.norm(cells - [hx, hy], axis=1) < 0.6
mix = np.bincount(ct[inside], minlength=3) / inside.sum()
x0 = 2.6; y0 = -0.75
a.text(x0 - 0.15, y0, "highlighted spot reads:", ha="right", va="center", fontsize=8.4, color=TEXT2)
left = x0
for k in range(3):
    a.add_patch(mp.Rectangle((left, y0 - 0.12), 2.2 * mix[k], 0.28, fc=TC[k], ec="white", lw=0.8, clip_on=False))
    left += 2.2 * mix[k]
a.text(x0 + 2.35, y0, f"a mix of {inside.sum()} cells", ha="left", va="center", fontsize=8.4, color=CRIMSON)
a.set_title("(a) Spot-based (Visium)", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=40)
a.text(0, 1.1, "whole transcriptome; each 55 µm spot covers several cells", transform=a.transAxes, fontsize=8.8, color=TEXT2)
# (b) imaging
b = axes[1]
gcol = [CRIMSON, "#2B2F36", tint(0.5)]
for (x, y), t in zip(cells, ct):
    n = rng.poisson(4)
    for _ in range(n):
        r, th = 0.14 * np.sqrt(rng.uniform()), rng.uniform(0, 2 * np.pi)
        g = t if rng.uniform() < 0.8 else rng.integers(0, 3)
        b.scatter([x + r * np.cos(th)], [y + r * np.sin(th)], s=3.5, color=gcol[[1, 0, 2][g]], lw=0, zorder=3)
b.set_title("(b) Imaging-based (MERFISH, Xenium, CosMx)", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=40)
b.text(0, 1.1, "each dot = one transcript, placed inside its cell; a fixed panel of genes", transform=b.transAxes,
       fontsize=8.8, color=TEXT2)
fig.text(0.125, 0.1, "Same trade-off as the rest of this chapter: resolution (cell or finer) vs breadth (all genes vs a "
         "panel of hundreds to a few thousand)", fontsize=8.9, color=TEXT2)
header(a, "Single-cell RNA-seq · spatial", "Dissociation throws away position; spatial methods keep it", dy=58, x1=2.1)
fig.savefig(OUT + "spatial.png"); print("ok k28 spatial", inside.sum())
