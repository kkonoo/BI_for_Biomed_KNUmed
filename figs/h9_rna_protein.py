from gallery import *
rng = np.random.default_rng(9)

# mRNA level only loosely predicts protein level across genes (r ≈ 0.4)
fig = plt.figure(figsize=(9.8, 4.9))
ax = fig.add_axes([0.07, 0.12, 0.5, 0.66])
DARK = "#5E6670"; GREYF = "#B8BEC7"
n = 900
x = rng.normal(0, 1, n)
y = 0.4 * x + np.sqrt(1 - 0.4 ** 2) * rng.normal(0, 1, n)
r = np.corrcoef(x, y)[0, 1]
mx, my = 1 + 0.7 * x, 3 + 0.9 * y                      # log10 mRNA (TPM), log10 protein (copies)
ax.scatter(mx, my, s=9, color=GREYF, lw=0, alpha=0.8)
b = np.polyfit(mx, my, 1); xs = np.linspace(mx.min(), mx.max(), 10)
ax.plot(xs, np.polyval(b, xs), color=TEXT2, lw=1.2, ls=(0, (4, 2)))
g1, g2 = (1.6, 4.85), (1.6, 3.35)                       # same mRNA, ~30× different protein
for (gx, gy), lab in ((g1, "gene P"), (g2, "gene Q")):
    ax.scatter([gx], [gy], s=46, color=CRIMSON, lw=1.2, edgecolor="white", zorder=3)
    ax.text(gx + 0.12, gy, lab, ha="left", va="center", fontsize=9.5, fontweight="semibold", color=CRIMSON)
ax.annotate("", (1.48, g2[1] + 0.05), (1.48, g1[1] - 0.05),
            arrowprops=dict(arrowstyle="<->", color=CRIMSON, lw=1))
ax.text(1.38, (g1[1] + g2[1]) / 2, "same mRNA,\n~30× protein", ha="right", va="center", fontsize=8.8,
        color=CRIMSON, linespacing=1.25, bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.5))
ax.text(0.02, 0.97, f"r ≈ {r:.1f}", transform=ax.transAxes, ha="left", va="top", fontsize=12, fontweight="bold",
        color=TEXT)
ax.set_xlabel("mRNA level (log scale)", color=TEXT2, fontsize=9.5)
ax.set_ylabel("protein level (log scale)", color=TEXT2, fontsize=9.5)
ax.set_xticks([]); ax.set_yticks([])
for sp in ("left", "bottom"): ax.spines[sp].set_color(LINE)
ax.set_title("One dot = one gene", loc="left", fontsize=12, fontweight="semibold", color=TEXT, pad=8)

# what sits between mRNA and protein
t = fig.add_axes([0.62, 0.12, 0.37, 0.66]); t.axis("off"); t.set_xlim(0, 1); t.set_ylim(0, 1)
t.text(0, 1.0, "What RNA-seq cannot see", ha="left", va="top", fontsize=12, fontweight="semibold", color=TEXT)
rows = [("translation rate", "how many proteins each mRNA makes"),
        ("protein half-life", "minutes for some, weeks for others"),
        ("post-translational control", "modification, localisation, cleavage")]
for k, (h, d) in enumerate(rows):
    y = 0.78 - k * 0.22
    t.scatter([0.02], [y], marker="s", s=22, color=tint(0.55), lw=0)
    t.text(0.07, y, h, ha="left", va="center", fontsize=10.3, fontweight="semibold", color=TEXT)
    t.text(0.07, y - 0.08, d, ha="left", va="center", fontsize=9, color=TEXT2)
t.text(0, 0.06, "→ measure protein directly (W6)", ha="left", va="center", fontsize=9.6, color=CRIMSON,
       fontweight="semibold")
header(ax, "Bulk RNA-seq · RNA is not protein", "mRNA level only loosely predicts protein level", dy=20, x1=1.84)
fig.savefig(OUT + "rna_protein.png"); print("ok h9 rna protein", round(r, 2))
