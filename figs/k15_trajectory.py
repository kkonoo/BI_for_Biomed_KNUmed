from gallery import *
from matplotlib.colors import LinearSegmentedColormap
rng = np.random.default_rng(15)

# A trajectory algorithm always returns a trajectory: real process vs no process vs a flipped root
fig, axes = plt.subplots(1, 3, figsize=(9.8, 3.9), gridspec_kw={"wspace": 0.08})
CM = LinearSegmentedColormap.from_list("pt", ["#D5D9DF", tint(0.45), CRIMSON])
def fit_curve(P, deg=3):
    """Principal-curve stand-in: cubic fit of y on x, then arc-length pseudotime of each cell's nearest curve point."""
    b = np.polyfit(P[:, 0], P[:, 1], deg); xs = np.linspace(P[:, 0].min(), P[:, 0].max(), 300)
    C = np.c_[xs, np.polyval(b, xs)]; s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(C, axis=0), axis=1))]
    pt = s[np.argmin(np.linalg.norm(P[:, None] - C[None], axis=2), axis=1)]
    return C, np.argsort(np.argsort(pt)) / (len(pt) - 1)          # rank-scaled, as tools report it
t = rng.uniform(0, 1, 260); real = np.c_[-2 + 4 * t, 1.2 * np.sin(2.4 * t - 0.6)] + rng.normal(0, 0.17, (260, 2))
blob = rng.normal(0, [0.85, 0.75], (260, 2))
cases = [(real, False, "✓ A real differentiation process", "pseudotime here means something", TEXT),
         (blob, False, "✕ Random cells, no process", "same algorithm, same confident output", CRIMSON),
         (real, True, "✕ Root chosen at the other end", "flip the start and the order flips", CRIMSON)]
for ax, (P, flip, title, sub, tc) in zip(axes, cases):
    C, pt = fit_curve(P, 3 if P is real else 2)
    if flip: pt, C = 1 - pt, C[::-1]
    ax.scatter(P[:, 0], P[:, 1], s=9, c=pt, cmap=CM, lw=0, alpha=0.9)
    ax.plot(C[:, 0], C[:, 1], color=TEXT, lw=1.4)
    ax.annotate("", C[-1], C[-12], arrowprops=dict(arrowstyle="-|>", color=TEXT, lw=1.4, mutation_scale=11))
    ax.scatter(*C[0], s=60, color=TEXT, zorder=4)
    ax.annotate("start (you choose it)" if flip else "start", C[0],
                xytext=(C[0][0] - 0.4, C[0][1] + 0.55) if flip else (C[0][0] + 0.1, C[0][1] - 0.75),
                fontsize=8.5, color=TEXT if not flip else CRIMSON, ha="center")
    ax.set_xlim(-2.9, 2.9); ax.set_ylim(-2.6, 2.4); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(True); sp.set_color("#E3E6EB")
    ax.set_title(title, loc="left", fontsize=11.8, fontweight="semibold", color=tc, pad=22)
    ax.text(0, 1.03, sub, transform=ax.transAxes, fontsize=8.9, color=TEXT2)
fig.text(0.125, -0.02, "colour = pseudotime (early grey → late red).  Pseudotime is an order, not time: it measures neither rate "
         "nor duration.  Test with a time course or lineage tracing.", fontsize=8.8, color=TEXT2)
header(axes[0], "Single-cell RNA-seq · trajectories", "A trajectory algorithm returns a trajectory, whether or not one exists",
       dy=40, x1=3.17)
fig.savefig(OUT + "trajectory_caveats.png"); print("ok k15 trajectory")
