from gallery import *

# Quantification: (a) when the samples are combined; (b) why TMT ratios are compressed; (c) by how much
fig = plt.figure(figsize=(9.8, 7.0))
gs = fig.add_gridspec(2, 2, height_ratios=[0.85, 1.05], width_ratios=[1, 1.15], hspace=0.48, wspace=0.3,
                      left=0.03, right=0.97, top=0.86, bottom=0.08)
DARK = "#5E6670"; GREYF = "#B8BEC7"

# (a) mixing point
a = fig.add_subplot(gs[0, :]); a.set_xlim(0, 10); a.set_ylim(0, 3.2); a.axis("off")
stages = ["cells", "protein", "peptides", "LC-MS run"]
SX = [2.6, 4.6, 6.6, 8.6]
for x, s in zip(SX, stages):
    a.text(x, 3.05, s, ha="center", va="center", fontsize=9.2, color=TEXT2)
    a.plot([x, x], [0.15, 2.85], color=FAINT, lw=8, zorder=0, solid_capstyle="butt")
rows = [("SILAC", 0, "heavy amino acids in culture"), ("TMT", 2.5, "chemical tag after digestion"),
        ("LFQ", None, "never mixed: one run each")]
for k, (nm, mix, sub) in enumerate(rows):
    y = 2.35 - k * 0.95
    a.text(0, y + 0.1, nm, ha="left", va="center", fontsize=11, fontweight="bold", color=TEXT)
    a.text(0, y - 0.18, sub, ha="left", va="center", fontsize=8.3, color=TEXT2)
    x0, x1 = 2.2, 9.2
    xm = SX[0] if mix == 0 else (SX[2] + 1.0 if mix == 2.5 else None)
    if xm is None:
        a.plot([x0, x1], [y + 0.1, y + 0.1], color=tint(0.6), lw=2.5, solid_capstyle="round")
        a.plot([x0, x1], [y - 0.1, y - 0.1], color=DARK, lw=2.5, solid_capstyle="round")
        a.text(x1 + 0.12, y, "every step\ncan drift", ha="left", va="center", fontsize=8.3, color=TEXT2,
               linespacing=1.15)
    else:
        if xm > x0:
            a.plot([x0, xm - 0.25], [y + 0.1, y + 0.1], color=tint(0.6), lw=2.5, solid_capstyle="round")
            a.plot([x0, xm - 0.25], [y - 0.1, y - 0.1], color=DARK, lw=2.5, solid_capstyle="round")
            a.plot([xm - 0.25, xm], [y + 0.1, y], color=tint(0.6), lw=2.5)
            a.plot([xm - 0.25, xm], [y - 0.1, y], color=DARK, lw=2.5)
        a.plot([xm, x1], [y, y], color=CRIMSON, lw=6, solid_capstyle="round")
        a.plot([xm, xm], [y + 0.28, y - 0.28], color=CRIMSON, lw=1)
        a.text(xm, y + 0.38, "mix", ha="center", va="bottom", fontsize=8.5, color=CRIMSON, fontweight="semibold")
        a.text(x1 + 0.12, y, "shared\nfrom here", ha="left", va="center", fontsize=8.3, color=CRIMSON,
               linespacing=1.15)
a.text(2.2, 0.05, "two samples (pink, grey) · thick red = both go through the same steps together", ha="left",
       va="center", fontsize=8.5, color=TEXT2)
a.set_title("(a) When are the samples combined?", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=8)

# (b) mechanism of ratio compression: reporter-ion intensities in two channels
b = fig.add_subplot(gs[1, 0])
t, bg, R = 1.0, 2.5, 8
X = np.array([0, 1])
b.bar(X, [R * t, t], color=CRIMSON, width=0.5, label="target peptide")
b.bar(X, [bg, bg], bottom=[R * t, t], color=GREYF, width=0.5, label="co-isolated background")
for x, tot, tg in zip(X, [R * t + bg, t + bg], [R * t, t]):
    b.text(x, tot + 0.25, f"{tot:g}", ha="center", va="bottom", fontsize=10, fontweight="bold", color=TEXT)
b.set_xticks(X); b.set_xticklabels(["channel 1", "channel 2"], fontsize=9, color=TEXT)
b.tick_params(axis="x", length=0); b.set_yticks([])
for sp in ("left", "top", "right"): b.spines[sp].set_visible(False)
b.spines["bottom"].set_color(LINE)
b.set_ylim(0, 13.5); b.set_xlim(-0.6, 2.6)
b.text(1.55, 9.5, f"true  {R * t:g} : {t:g} = {R}×", ha="left", va="center", fontsize=9.5, color=CRIMSON,
       fontweight="semibold")
b.text(1.55, 8.1, f"seen  {R * t + bg:g} : {t + bg:g} = {(R * t + bg) / (t + bg):.0f}×", ha="left", va="center",
       fontsize=9.5, color=TEXT, fontweight="semibold")
b.legend(frameon=False, fontsize=8.4, loc="upper right", bbox_to_anchor=(1.02, 0.55), labelcolor=TEXT,
         handlelength=1)
b.set_title("(b) Why TMT ratios shrink", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
b.text(0, 1.03, "background peptides add the same signal to every channel", transform=b.transAxes, fontsize=9,
       color=TEXT2)

# (c) true vs measured ratio
c = fig.add_subplot(gs[1, 1])
tr = np.logspace(0, np.log10(16), 100)
for bgr, col, lab, lw in ((0, GREYF, "0  (ideal)", 1.4), (0.25, DARK, "0.25  (e.g. SPS-MS3)", 1.8),
                          (2.5, CRIMSON, "2.5  (MS2, complex sample)", 2.4)):
    c.plot(tr, (tr + bgr) / (1 + bgr), color=col, lw=lw, label=lab, ls=(0, (4, 2)) if bgr == 0 else "-")
c.set_xscale("log"); c.set_yscale("log")
ticks = [1, 2, 4, 8, 16]
c.set_xticks(ticks); c.set_xticklabels([str(v) for v in ticks]); c.set_yticks(ticks); c.set_yticklabels([str(v) for v in ticks])
c.minorticks_off()
c.plot([8, 8], [1, 3], color=CRIMSON, lw=0.8, ls=":"); c.plot([1, 8], [3, 3], color=CRIMSON, lw=0.8, ls=":")
c.text(8.4, 2.6, "8× reads as 3×", ha="left", va="center", fontsize=9, color=CRIMSON, fontweight="semibold")
c.set_xlabel("true ratio", color=TEXT2, fontsize=9); c.set_ylabel("measured ratio", color=TEXT2, fontsize=9)
for sp in ("left", "bottom"): c.spines[sp].set_color(LINE)
c.tick_params(colors=TEXT2, labelsize=8.3)
c.legend(frameon=False, fontsize=8.4, loc="upper left", labelcolor=TEXT, title="background : target",
         title_fontsize=8.4, alignment="left")
c.set_title("(c) The bigger the change, the more it shrinks", loc="left", fontsize=12.3, fontweight="semibold",
            color=TEXT, pad=22)
c.text(0, 1.03, "measured = (true + bg) / (1 + bg)", transform=c.transAxes, fontsize=9, color=TEXT2)
header(a, "Proteomics · quantification", "Mix early to share the noise; TMT pays with compressed ratios",
       dy=26)
fig.savefig(OUT + "quantification.png"); print("ok m5 quant")
