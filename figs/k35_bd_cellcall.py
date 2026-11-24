from gallery import *
from scipy.ndimage import gaussian_filter1d
rng = np.random.default_rng(9)

# BD Rhapsody's basic cell calling: rank cell labels by reads, draw the cumulative reads curve on a log rank axis,
# and cut where its second derivative is lowest (the sharpest bend). Same simulated barcodes as k9_knee.py.
big = np.exp(rng.normal(np.log(6000), 0.55, 6000))
small = np.exp(rng.normal(np.log(420), 0.35, 1200))
empty = np.concatenate([np.exp(rng.normal(np.log(30), 0.75, 150000)), np.exp(rng.normal(np.log(260), 0.4, 2500))])
reads = np.concatenate([big, small, empty]); cell = np.r_[np.ones(big.size + small.size), np.zeros(empty.size)]
o = np.argsort(-reads); reads, cell = reads[o], cell[o]
rank = np.arange(1, reads.size + 1); cum = np.cumsum(reads) / reads.sum()
g = np.linspace(0, np.log10(rank[-1]), 400)
c = gaussian_filter1d(np.interp(g, np.log10(rank), cum), 4)
d2 = np.gradient(np.gradient(c, g), g)
i = int(np.argmin(d2)); r_cut = 10 ** g[i]
missed = int(((cell == 1) & (rank > r_cut)).sum())
DARK = "#5E6670"; GREYF = "#B8BEC7"
fig, axes = plt.subplots(1, 2, figsize=(9.8, 3.7), gridspec_kw={"wspace": 0.28})
fig.subplots_adjust(left=0.07, right=0.98, top=0.74, bottom=0.17)
a, b = axes
a.plot(rank, cum, color=GREYF, lw=2.2)
ci = cell == 1
a.plot(rank[ci & (rank <= r_cut)], cum[ci & (rank <= r_cut)], color=CRIMSON, lw=2.6)
for ax in axes:
    ax.axvline(r_cut, color=TEXT, lw=0.9, ls=(0, (4, 3)))
    ax.set_xscale("log"); ax.set_xlim(1, 3e5); ax.minorticks_off()
    for sp in ("left", "bottom"): ax.spines[sp].set_color(LINE)
    ax.tick_params(colors=TEXT2, labelsize=8.3)
    ax.set_xlabel("cell labels, ranked by reads", color=TEXT2, fontsize=9)
a.set_ylim(0, 1.02); a.set_yticks([0, 0.5, 1]); a.set_yticklabels(["0", "50%", "100%"])
a.set_ylabel("cumulative share of reads", color=TEXT2, fontsize=9)
a.text(r_cut * 1.2, 0.3, f"cut at rank {r_cut:,.0f}\n→ putative cells on the left", ha="left", va="center",
       fontsize=8.6, color=TEXT, linespacing=1.25)
a.text(3, 0.62, "rises steeply:\neach label adds\nmany reads", ha="left", va="center", fontsize=8.3, color=CRIMSON,
       linespacing=1.2)
a.set_title("(a) Cumulative reads curve", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
a.text(0, 1.03, "steep while labels are cells, flat once they are empty wells", transform=a.transAxes, fontsize=9,
       color=TEXT2)
b.plot(10 ** g, d2, color=DARK, lw=1.6)
b.axhline(0, color=LINE, lw=0.8)
b.scatter([r_cut], [d2[i]], s=40, color=CRIMSON, zorder=3)
b.annotate("minimum = the sharpest bend", (r_cut, d2[i]), xytext=(25, d2[i] * 0.8), fontsize=8.6,
           color=CRIMSON, fontweight="semibold", arrowprops=dict(arrowstyle="-", color=CRIMSON, lw=0.7))
b.set_yticks([0]); b.set_ylabel("second derivative", color=TEXT2, fontsize=9)
b.set_title("(b) Its second derivative", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
b.text(0, 1.03, f"lowest point = the cut · {missed:,} small real cells fall right of it", transform=b.transAxes,
       fontsize=9, color=TEXT2)
header(a, "Single-cell RNA-seq · BD Rhapsody cell calling", "Cut where the cumulative reads curve bends most",
       dy=40, x1=2.3)
fig.savefig(OUT + "bd_cell_calling.png"); print("ok k35 bd", round(r_cut), missed)
