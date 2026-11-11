from gallery import *
rng = np.random.default_rng(16)

# Phosphoproteomics: (a) phosphopeptides are rare, so they are enriched first; (b) more phosphopeptide can just be
# more protein; (c) which residue carries the phosphate is itself an estimate
fig, axes = plt.subplots(1, 3, figsize=(9.8, 4.2), gridspec_kw={"wspace": 0.3, "width_ratios": [1.05, 0.9, 1.05]})
DARK = "#5E6670"; GREYF = "#B8BEC7"

def title(ax, t, sub):
    ax.set_title(t, loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
    ax.text(0, 1.03, sub, transform=ax.transAxes, fontsize=9, color=TEXT2)

# (a) enrichment: 100 peptides before and after (illustrative shares)
a = axes[0]; a.set_xlim(0, 10); a.set_ylim(0, 10); a.axis("off"); a.set_aspect("equal", adjustable="datalim")
for k, (nph, lab) in enumerate(((2, "digest"), (85, "after TiO₂ / IMAC"))):
    x0 = 0.3 + k * 5.0
    ph = set(rng.choice(100, nph, replace=False))
    for i in range(100):
        a.add_patch(mp.Circle((x0 + (i % 10) * 0.42 + 0.2, 7.6 - (i // 10) * 0.42), 0.16,
                              fc=CRIMSON if i in ph else GREYF, ec="none"))
    a.text(x0 + 2.1, 2.9, lab, ha="center", va="center", fontsize=9, color=TEXT)
    a.text(x0 + 2.1, 2.2, f"~{nph}% phospho", ha="center", va="center", fontsize=9.5, fontweight="bold",
           color=CRIMSON if k else TEXT2)
a.annotate("", (5.1, 5.7), (4.6, 5.7), arrowprops=dict(arrowstyle="-|>", color=TEXT2, lw=1, mutation_scale=8))
title(a, "(a) Enrich first", "red = phosphopeptide (illustrative)")

# (b) occupancy: both phosphopeptide and protein double
b = axes[1]
vals = [("phospho-\npeptide", 1, 2), ("whole\nprotein", 5, 10), ("phospho ÷\nprotein", 0.2, 0.2)]
X = np.arange(3)
for x, (nm, v1, v2) in zip(X, vals):                      # bars scaled so that the treated value = 1
    b.bar(x - 0.17, v1 / v2, width=0.32, color=GREYF)
    b.bar(x + 0.17, 1, width=0.32, color=CRIMSON)
    b.text(x + 0.17, 1.04, f"×{v2 / v1:g}", ha="center", va="bottom", fontsize=9.5, fontweight="bold",
           color=CRIMSON if v2 != v1 else TEXT)
b.set_xticks(X); b.set_xticklabels([v[0] for v in vals], fontsize=8.6, color=TEXT)
b.tick_params(axis="x", length=0); b.set_yticks([]); b.set_ylim(0, 1.3)
for sp in ("left", "top", "right"): b.spines[sp].set_visible(False)
b.spines["bottom"].set_color(LINE)
b.text(0.0, 1.0, "grey = control\nred = treated", transform=b.transAxes, ha="left", va="top", fontsize=8.3,
       color=TEXT2, linespacing=1.2)
title(b, "(b) Or just more protein?", "occupancy did not change")

# (c) site localisation
c = axes[2]; c.set_xlim(0, 10); c.set_ylim(0, 10); c.axis("off")
pep = list("AGSPTSLAK")
cand = {2: 0.62, 4: 0.30, 5: 0.08}
for i, aa in enumerate(pep):
    x = 0.6 + i * 1.0
    on = i in cand
    c.add_patch(mp.FancyBboxPatch((x - 0.38, 7.3), 0.76, 0.9, boxstyle="round,pad=0,rounding_size=0.12",
                                  fc=DARK if on else FAINT, ec="none"))
    c.text(x, 7.75, aa, ha="center", va="center", fontsize=10.5, fontweight="bold", color="white" if on else TEXT2)
    if on:
        p = cand[i]
        c.add_patch(mp.Rectangle((x - 0.3, 2.0), 0.6, 4.4 * p, fc=CRIMSON if p >= 0.75 else tint(0.55), ec="none"))
        c.text(x, 2.0 + 4.4 * p + 0.15, f"{p:.2f}", ha="center", va="bottom", fontsize=9, fontweight="bold", color=TEXT)
c.plot([0.1, 9.4], [2.0 + 4.4 * 0.75] * 2, color=TEXT, lw=0.8, ls=(0, (3, 2)))
c.text(9.4, 2.0 + 4.4 * 0.75 + 0.2, "0.75", ha="right", va="bottom", fontsize=8.3, color=TEXT)
c.plot([0.1, 9.4], [2.0, 2.0], color=LINE, lw=0.8)
c.text(0.1, 1.3, "probability that each S/T/Y carries the phosphate:\nnone reaches 0.75 → site uncertain", ha="left",
       va="center", fontsize=8.6, color=TEXT2, linespacing=1.2)
title(c, "(c) Which residue?", "one phosphopeptide, three candidate sites")
header(axes[0], "Proteomics · phosphorylation", "A phosphosite is enriched, normalised, and localised — each step can mislead",
       dy=40, x1=3.75)
fig.savefig(OUT + "phospho.png"); print("ok m16 phospho")
