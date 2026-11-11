from gallery import *

# Three things besides biology that change a raw count: depth, gene length, composition
fig = plt.figure(figsize=(9.8, 6.0))
gs = fig.add_gridspec(2, 3, height_ratios=[1, 0.32], hspace=0.55, wspace=0.28,
                      left=0.06, right=0.99, top=0.8, bottom=0.03)
DARK = "#5E6670"; GREYF = "#B8BEC7"

def panel(ax, t, sub):
    ax.set_title(t, loc="left", fontsize=12, fontweight="semibold", color=TEXT, pad=22)
    ax.text(0, 1.03, sub, transform=ax.transAxes, ha="left", va="bottom", fontsize=9, color=TEXT2)
def style(ax):
    for sp in ("left", "bottom"): ax.spines[sp].set_color(LINE)
    ax.tick_params(colors=TEXT2, labelsize=8.5)

# ① depth
a = fig.add_subplot(gs[0, 0]); style(a)
panel(a, "① Library size (depth)", "same gene, same expression")
a.bar([0, 1], [100, 200], width=0.55, color=[GREYF, tint(0.55)])
for x, v in ((0, 100), (1, 200)): a.text(x, v + 6, f"{v}", ha="center", va="bottom", fontsize=10, fontweight="semibold", color=TEXT)
a.set_xticks([0, 1]); a.set_xticklabels(["sample A\n10 M reads", "sample B\n20 M reads"], fontsize=8.8)
a.set_ylim(0, 260); a.set_ylabel("raw count", color=TEXT2, fontsize=9)
a.text(0.5, 235, "2× deeper → 2× the count", ha="center", va="center", fontsize=8.8, color=CRIMSON, fontweight="semibold")

# ② gene length
b = fig.add_subplot(gs[0, 1]); b.set_xlim(0, 10); b.set_ylim(0, 10); b.axis("off")
panel(b, "② Gene length", "same number of molecules")
def molecules(y0, L, lab, col):
    b.text(0, y0 + 1.75, lab, ha="left", va="center", fontsize=9, fontweight="semibold", color=TEXT)
    for k in range(3):
        y = y0 + 1.0 - k * 0.55
        b.plot([0.2, 0.2 + L], [y, y], color=LINE, lw=3, solid_capstyle="butt")
        cuts = np.arange(0.2, 0.2 + L - 0.05, 0.95)
        for c in cuts:
            b.plot([c + 0.06, min(c + 0.86, 0.2 + L - 0.02)], [y, y], color=col, lw=3, solid_capstyle="butt")
molecules(6.3, 1.9, "short transcript (1 kb) × 3", GREYF)
molecules(1.9, 7.6, "long transcript (4 kb) × 3", tint(0.55))
b.text(2.6, 7.3, "→ 6 fragments", ha="left", va="center", fontsize=9, color=TEXT2)
b.text(0, 0.3, "→ 24 fragments: 4× the reads", ha="left", va="center", fontsize=9, color=CRIMSON, fontweight="semibold")

# ③ composition
c = fig.add_subplot(gs[0, 2]); c.set_xlim(0, 10); c.set_ylim(0, 10); c.axis("off")
panel(c, "③ Composition", "same depth, one gene takes over")
def library(y, shares, lab):
    c.text(0, y + 1.0, lab, ha="left", va="center", fontsize=9, fontweight="semibold", color=TEXT)
    x = 0
    for k, (sh, col, t) in enumerate(shares):
        c.add_patch(mp.Rectangle((x, y - 0.45), 10 * sh, 0.9, fc=col, ec="white", lw=1.2))
        if t: c.text(x + 5 * sh, y, t, ha="center", va="center", fontsize=8.5, color="white" if col in (CRIMSON, DARK) else TEXT,
                     fontweight="semibold")
        x += 10 * sh
genes = [("A", "#D5D9DF"), ("B", "#C9CDD3"), ("C", "#D5D9DF"), ("D", "#C9CDD3")]
library(7.0, [(0.05, CRIMSON, "")] + [(0.2375, col, g) for g, col in genes], "control")
library(3.6, [(0.6, CRIMSON, "gene X")] + [(0.1, col, g) for g, col in genes], "treated")
c.text(0, 1.85, "gene A: 24% → 10% of the library", ha="left", va="center", fontsize=8.8, color=TEXT)
c.text(0, 0.95, "→ looks 2.4× down, but did not change", ha="left", va="center", fontsize=8.8, color=CRIMSON,
       fontweight="semibold")

# what fixes what
t = fig.add_subplot(gs[1, :]); t.set_xlim(0, 10); t.set_ylim(0, 1); t.axis("off")
units = [("Raw count", "nothing", "input to DESeq2 / edgeR", False),
         ("CPM", "①", "filter low counts", False),
         ("TPM", "① + ②", "report, heatmaps", False),
         ("TMM / RLE", "① + ③", "the DE standard", True)]
for k, (u, fix, use, hl) in enumerate(units):
    x = k * 2.52
    t.add_patch(mp.FancyBboxPatch((x, 0.05), 2.38, 0.9, boxstyle="round,pad=0,rounding_size=0.06",
                                  fc=tint(0.1) if hl else "#F8F9FA", ec=CRIMSON if hl else LINE, lw=1 if hl else 0.9))
    t.text(x + 0.15, 0.68, u, ha="left", va="center", fontsize=10.5, fontweight="bold", color=CRIMSON if hl else TEXT)
    t.text(x + 2.25, 0.68, "fixes " + fix, ha="right", va="center", fontsize=8.8, color=TEXT2)
    t.text(x + 0.15, 0.3, use, ha="left", va="center", fontsize=8.8, color=TEXT)

header(a, "Bulk RNA-seq · normalisation", "A raw count is not just biology", dy=40, x1=3.6)
fig.savefig(OUT + "normalization.png"); print("ok h2 normalisation")
