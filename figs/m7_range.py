from gallery import *

# Plasma dynamic range: typical concentrations in healthy adults (order of magnitude), what each method covers,
# and how the protein mass is shared out
fig = plt.figure(figsize=(9.8, 5.0))
gs = fig.add_gridspec(1, 2, width_ratios=[1, 0.24], wspace=0.25, left=0.2, right=0.93, top=0.8, bottom=0.12)
DARK = "#5E6670"; GREYF = "#B8BEC7"
prot = [("albumin", 10.6), ("immunoglobulin G", 10.0), ("fibrinogen", 9.5), ("transferrin", 9.4),
        ("C-reactive protein", 6.0), ("PSA", 3.0), ("troponin I", 0.5), ("IL-6", 0.3), ("TNF-α", 0.0)]
a = fig.add_subplot(gs[0])
Y = np.arange(len(prot))
for y, (nm, v) in zip(Y, prot):
    top = v > 9
    a.plot([-1, v], [y, y], color=DARK if top else GREYF, lw=7, solid_capstyle="butt")
    a.text(-1.2, y, nm, ha="right", va="center", fontsize=9.5, color=TEXT)
a.invert_yaxis()
a.set_xlim(-1, 11.2); a.set_ylim(len(prot) + 1.3, -1.6)
a.set_xticks(range(0, 12, 2))
a.set_xticklabels(["1 pg", "100 pg", "10 ng", "1 µg", "100 µg", "10 mg"], fontsize=8.5, color=TEXT2)
a.set_yticks([])
for sp in ("left",): a.spines[sp].set_visible(False)
a.spines["bottom"].set_color(LINE); a.tick_params(axis="x", colors=TEXT2)
a.set_xlabel("concentration per mL of plasma (log scale)", color=TEXT2, fontsize=9)
# coverage windows
ymS, ymA = -0.95, len(prot) + 0.6
a.annotate("", (10.9, ymS), (5.2, ymS), arrowprops=dict(arrowstyle="<->", color=DARK, lw=1.2))
a.text(5.0, ymS, "MS without depletion: ~5–6 orders", ha="right", va="center", fontsize=8.8, color=DARK,
       fontweight="semibold")
a.annotate("", (-0.8, ymA), (10.9, ymA), arrowprops=dict(arrowstyle="->", color=CRIMSON, lw=1.4))
a.text(5.0, ymA - 0.5, "affinity panels (Olink, SomaScan) reach down to pg/mL", ha="center", va="center",
       fontsize=8.8, color=CRIMSON, fontweight="semibold")
a.annotate("", (0.3, 7.5), (10.6, 7.5), arrowprops=dict(arrowstyle="<->", color=TEXT2, lw=0.8))
a.text(5.45, 7.5, "> 10 orders of magnitude", ha="center", va="center", fontsize=8.8, color=TEXT2,
       bbox=dict(fc="white", ec="none", pad=1.5))
a.set_title("Typical levels in healthy adults", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22,
            x=-0.22)
a.text(-0.22, 1.03, "order of magnitude only; biomarkers sit at the bottom", transform=a.transAxes, fontsize=9,
       color=TEXT2)

# share of protein mass
b = fig.add_subplot(gs[1])
parts = [("albumin\n≈ half", 55, DARK), ("immuno-\nglobulins", 18, "#8A929C"), ("18 more\nabundant", 26, GREYF),
         ("all the\nothers ≈ 1%", 1, CRIMSON)]
bottom = 0
for nm, v, c in parts:
    b.bar(0, v, bottom=bottom, color=c, width=0.6, edgecolor="white", lw=1)
    yy = bottom + v / 2 if v > 5 else bottom + v + 3
    b.text(0.4, yy, nm, ha="left", va="center", fontsize=8.6, color=CRIMSON if c == CRIMSON else TEXT,
           fontweight="semibold" if c == CRIMSON else "normal", linespacing=1.1)
    bottom += v
b.set_xlim(-0.4, 1.4); b.set_ylim(0, 106); b.axis("off")
b.set_title("Plasma protein mass", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
b.text(0, 1.03, "top 20 proteins ≈ 99%", transform=b.transAxes, fontsize=9, color=TEXT2)
header(a, "Proteomics · plasma", "The biomarkers you want sit ten orders below albumin", x=-0.22, dy=40, x1=1.5)
fig.savefig(OUT + "dynamic_range.png"); print("ok m7 range")
