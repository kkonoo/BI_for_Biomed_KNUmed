from gallery import *

# Droplet loading: cells are loaded dilute so most droplets hold a bead but no cell (Poisson);
# loading more cells per run raises the share of droplets with two cells
fig = plt.figure(figsize=(9.8, 4.4))
gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1], wspace=0.2, left=0.01, right=0.99, top=0.8, bottom=0.14)
DARK = "#5E6670"; GREYF = "#C3C8CF"
a = fig.add_subplot(gs[0]); a.set_xlim(0, 10); a.set_ylim(0, 6); a.axis("off")
a.set_title("(a) One droplet = one bead (+ maybe a cell)", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
a.text(0, 1.03, "every bead carries its own barcode", transform=a.transAxes, fontsize=9, color=TEXT2)
# inlets
for y, lab in ((5.1, "cells"), (3.0, "barcoded gel beads"), (0.9, "oil")):
    a.plot([0.2, 2.6], [y, y], color=LINE, lw=6, solid_capstyle="round", zorder=1)
    a.text(0.2, y + 0.35, lab, ha="left", va="center", fontsize=8.6, color=TEXT2)
for x in (0.6, 1.7): a.add_patch(mp.Circle((x, 5.1), 0.14, fc=DARK, ec="none", zorder=2))
for x in (0.5, 1.0, 1.5, 2.0): a.add_patch(mp.Circle((x, 3.0), 0.12, fc=tint(0.45), ec="none", zorder=2))
a.plot([2.6, 3.2], [5.1, 3.15], color=LINE, lw=6, solid_capstyle="round", zorder=1)
a.plot([2.6, 3.2], [0.9, 2.85], color=LINE, lw=6, solid_capstyle="round", zorder=1)
a.plot([2.6, 9.8], [3.0, 3.0], color=FAINT, lw=40, solid_capstyle="butt", zorder=0)
seq = [0, 0, 1, 0, 0, 0, 2, 0, 1, 0]
for i, n in enumerate(seq):
    x = 3.6 + i * 0.62
    a.add_patch(mp.Circle((x, 3.0), 0.27, fc="white", ec=CRIMSON if n == 2 else LINE, lw=1.2 if n == 2 else 0.9, zorder=2))
    a.add_patch(mp.Circle((x + 0.08, 2.9), 0.09, fc=tint(0.45), ec="none", zorder=3))
    for k in range(n): a.add_patch(mp.Circle((x - 0.08 + 0.12 * k, 3.09 - 0.05 * k), 0.09, fc=DARK, ec="none", zorder=3))
a.text(3.6 + 6 * 0.62, 3.55, "2 cells", ha="center", va="center", fontsize=8.3, color=CRIMSON, fontweight="semibold")
a.text(3.6 + 2 * 0.62, 2.45, "1 cell", ha="center", va="center", fontsize=8.3, color=TEXT)
a.text(3.6, 2.45, "empty", ha="center", va="center", fontsize=8.3, color=TEXT2)
a.text(3.3, 1.4, "Cells are loaded dilute, so that two cells rarely share a droplet.\nThe price: most droplets hold a bead and no cell —\nthey still get a barcode and still produce reads (§3.3).",
       ha="left", va="center", fontsize=8.6, color=TEXT2, linespacing=1.3)

# (b) Poisson loading
b = fig.add_subplot(gs[1])
lam = np.array([0.05, 0.1, 0.2, 0.4])
p0 = np.exp(-lam); p1 = lam * np.exp(-lam); occ = 1 - p0
multi = (occ - p1) / occ
X = np.arange(len(lam))
b.bar(X, 1 - multi, color=DARK, width=0.6, label="one cell")
b.bar(X, multi, bottom=1 - multi, color=CRIMSON, width=0.6, label="two or more cells")
for x, m, e in zip(X, multi, p0):
    b.text(x, 1.02, f"{100 * m:.0f}%", ha="center", va="bottom", fontsize=9.5, color=CRIMSON, fontweight="bold")
    b.text(x, -0.14, f"{100 * e:.0f}% empty", ha="center", va="top", fontsize=8.2, color=TEXT2)
b.set_xticks(X); b.set_xticklabels([f"λ = {v:g}" for v in lam], fontsize=8.6, color=TEXT)
b.tick_params(axis="x", length=0, pad=2)
b.set_ylim(0, 1.15); b.set_yticks([0, 0.5, 1]); b.set_yticklabels(["0", "50%", "100%"])
for sp in ("left", "bottom"): b.spines[sp].set_color(LINE)
b.tick_params(axis="y", colors=TEXT2, labelsize=8.3)
b.set_ylabel("droplets that hold a cell", color=TEXT2, fontsize=9)
b.legend(frameon=False, fontsize=8.4, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2, labelcolor=TEXT,
         handlelength=1)
b.set_title("(b) Load more cells, get more doublets", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
b.text(0, 1.03, "Poisson: λ = mean cells per droplet", transform=b.transAxes, fontsize=9, color=TEXT2)
header(a, "Single-cell RNA-seq · droplets", "Most droplets are empty on purpose", dy=40, x1=1.98)
fig.savefig(OUT + "droplet_loading.png"); print("ok k18 droplet", np.round(multi, 3))
