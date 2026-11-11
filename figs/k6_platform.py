from gallery import *

# The single-cell trade-off: number of cells vs depth per cell, on one sequencing budget line
fig, ax = plt.subplots(figsize=(9.8, 5.0))
DARK = "#5E6670"
pts = [("Bulk RNA-seq", 1, 3e7, DARK, "one averaged 'cell',\nvery deep", "above"),
       ("Plate\n(Smart-seq2/3)", 400, 3.2e6, tint(0.55), "full-length, isoforms;\nlabour-intensive", "diag"),
       ("Droplet\n(10x Genomics)", 50000, 26000, CRIMSON, "3′ or 5′ end only;\nthe current default", "diag"),
       ("Combinatorial\nindexing (sci-)", 1_000_000, 1300, tint(0.35), "whole-organism atlases;\nvery shallow",
        "diag")]
x = np.logspace(1.9, 6.75, 80)
ax.plot(x, 1.3e9 / x, color=LINE, lw=1.6, ls=(0, (5, 3)), zorder=1)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(0.3, 4e7); ax.set_ylim(3e2, 5e8)
fig.canvas.draw()
p1, p2 = ax.transData.transform([(1e3, 1.3e6), (1e4, 1.3e5)])
ang = np.degrees(np.arctan2(p2[1] - p1[1], p2[0] - p1[0]))
ax.text(4e3, 1.3e9 / 4e3 * 2.1, "one fixed sequencing budget", fontsize=9, color=TEXT2, rotation=ang,
        rotation_mode="anchor", ha="center", va="center")
for nm, n, depth, c, note, pos in pts:
    ax.scatter([n], [depth], s=330, color=c, alpha=0.18, lw=0, zorder=3)
    ax.scatter([n], [depth], s=80, color=c, lw=0, zorder=4)
    kw = dict(textcoords="offset points", fontsize=10.2, fontweight="bold", linespacing=1.15,
              color=CRIMSON if c == CRIMSON else TEXT)
    nk = dict(textcoords="offset points", fontsize=8.8, color=TEXT2, linespacing=1.2)
    if pos == "above":
        ax.annotate(nm, (n, depth), xytext=(0, 16), ha="center", va="bottom", **kw)
        ax.annotate(note, (n, depth), xytext=(0, -14), ha="center", va="top", **nk)
    else:
        ax.annotate(nm, (n, depth), xytext=(11, 7), ha="left", va="bottom", **kw)
        ax.annotate(note, (n, depth), xytext=(-11, -7), ha="right", va="top", **nk)
ax.set_xlabel("cells profiled per experiment", color=TEXT2, fontsize=9.5)
ax.set_ylabel("reads (or UMIs) per cell", color=TEXT2, fontsize=9.5)
for sp in ("left", "bottom"): ax.spines[sp].set_color(LINE)
ax.tick_params(colors=TEXT2, labelsize=8.5)
ax.minorticks_off()
header(ax, "Single-cell RNA-seq · platforms", "Every platform buys cells with depth, or depth with cells")
fig.savefig(OUT + "platform_tradeoff.png"); print("ok k6 platform")
