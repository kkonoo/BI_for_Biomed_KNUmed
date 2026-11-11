from gallery import *

# The kit-ome. Model: every extraction adds the same small amount of bacterial DNA from the reagents
# (K ≈ 10³ cell-equivalents, the order of magnitude seen in dilution experiments). Kit share = K / (K + sample).
# (a) kit share vs input biomass, with typical sample types; (b) one mock community diluted step by step
fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.5), gridspec_kw={"wspace": 0.3, "width_ratios": [1.15, 1]})
DARK = "#5E6670"; GREYF = "#B8BEC7"
K = 1e3
share = lambda cells: K / (K + cells)

def style(ax, t, sub):
    ax.set_title(t, loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
    ax.text(0, 1.03, sub, transform=ax.transAxes, fontsize=9, color=TEXT2)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"): ax.spines[sp].set_color(LINE)
    ax.tick_params(colors=TEXT2, labelsize=8.6)

a = axes[0]
x = np.logspace(0, 11, 300)
a.axhspan(50, 100, color=FAINT, lw=0, zorder=0)
a.plot(x, 100 * share(x), color=CRIMSON, lw=2.2, zorder=3)
a.text(10 ** 10.9, 93, "most reads come\nfrom the kit", ha="right", va="top", fontsize=8.6, color=TEXT2, linespacing=1.15)
a.set_xscale("log"); a.set_xlim(1, 1e11); a.set_ylim(-34, 104)
a.set_yticks([0, 50, 100]); a.set_yticklabels(["0", "50%", "100%"])
a.spines["left"].set_bounds(0, 100)
a.set_xticks([1e1, 1e3, 1e5, 1e7, 1e9, 1e11])
a.set_xlabel("bacterial cells in the sample", color=TEXT2, fontsize=9)
bands = [("blood, placenta,\ntumour tissue", 10 ** 0.5, 1e3, CRIMSON), ("skin swab,\nlung lavage", 1e3, 1e6, tint(0.5)),
         ("saliva,\nvaginal swab", 1e7, 1e9, GREYF), ("stool", 10 ** 9.5, 1e11, DARK)]
for k, (nm, lo, hi, c) in enumerate(bands):
    y = -12 if k % 2 == 0 else -24
    a.plot([lo, hi], [y, y], color=c, lw=5, solid_capstyle="butt")
lab =[("blood, placenta, tumour", -7.5), ("skin swab, lung lavage", -19.5), ("saliva, vaginal swab", -7.5),
       ("stool", -19.5)]
for (nm, lo, hi, c), (txt, ty) in zip(bands, lab):
    a.text(np.sqrt(lo * hi), ty, txt, ha="center", va="bottom", fontsize=7.8, color=TEXT)
a.text(1.2, -31, "typical input per extraction (rough orders of magnitude)", ha="left", va="center", fontsize=7.8,
       color=TEXT2)
a.set_ylabel("reads from kit DNA", color=TEXT2, fontsize=9)
a.yaxis.set_label_coords(-0.1, 0.62)
style(a, "(a) Less sample, more kit", "model: each extraction adds ~10³ cells' worth of DNA")

b = axes[1]
inputs = [1e8, 1e6, 1e5, 1e4, 1e3, 1e2]
kit = [("Ralstonia", 0.45, CRIMSON), ("Bradyrhizobium", 0.3, tint(0.6)), ("other kit genera", 0.25, tint(0.3))]
for i, c_ in enumerate(inputs):
    s = share(c_)
    b.bar(i, 1 - s, color=GREYF, width=0.66, edgecolor="white", lw=0.8)
    bot = 1 - s
    for nm, f, col in kit:
        b.bar(i, s * f, bottom=bot, color=col, width=0.66, edgecolor="white", lw=0.8)
        bot += s * f
    b.text(i, 1.02, f"{100 * s:.0f}%" if s >= 0.01 else "<1%", ha="center", va="bottom", fontsize=8.6,
           color=CRIMSON if s >= 0.5 else TEXT, fontweight="semibold")
b.set_xticks(range(6)); b.set_xticklabels(["10⁸", "10⁶", "10⁵", "10⁴", "10³", "10²"])
b.set_xlabel("cells of the mock community put in", color=TEXT2, fontsize=9)
b.set_ylim(0, 1.12); b.set_yticks([0, 0.5, 1]); b.set_yticklabels(["0", "50%", "100%"])
b.spines["left"].set_bounds(0, 1)
handles = [mp.Patch(fc=GREYF, label="the mock community")] + [mp.Patch(fc=col, label=nm) for nm, f, col in kit]
b.legend(handles=handles, frameon=False, fontsize=8.2, loc="upper left", bbox_to_anchor=(-0.02, -0.17), ncol=2,
         labelcolor=TEXT, handlelength=1, columnspacing=1)
style(b, "(b) Dilute a known community", "number on top = share of reads from the kit")
header(a, "Metagenomics · contamination", "With little biomass, most of what you sequence came out of the kit",
       dy=40, x1=2.2)
fig.savefig(OUT + "contamination.png"); print("ok n2")
