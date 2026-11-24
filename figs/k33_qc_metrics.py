from gallery import *
rng = np.random.default_rng(33)

# The first QC report (Cell Ranger web summary): where the reads go (valid barcode → on a gene → in a cell),
# how saturated the library is (reads that only re-read molecules already seen), and genes per cell
fig = plt.figure(figsize=(9.8, 3.9))
gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 1, 1], wspace=0.3, left=0.01, right=0.985, top=0.76, bottom=0.17)
DARK = "#5E6670"; GREYF = "#C3C8CF"

# (a) the read funnel
a = fig.add_subplot(gs[0]); a.set_xlim(0, 10); a.set_ylim(0, 4.6); a.axis("off")
a.set_title("(a) Where the reads go", loc="left", fontsize=12, fontweight="semibold", color=TEXT, pad=22)
a.text(0, 1.03, "share of all reads that survive each step", transform=a.transAxes, fontsize=8.8, color=TEXT2)
steps = [("all reads", None, ""),
         ("valid barcode", 0.96, "barcode not on the list"),
         ("mapped to a gene", 0.80, "intron, intergenic, unmapped"),
         ("in a cell", 0.88, "empty droplets (soup)")]
keep, L = 1.0, 6.0
for i, (nm, f, lost) in enumerate(steps):
    y = 3.95 - i * 1.05
    if f: keep *= f
    a.add_patch(mp.Rectangle((0, y - 0.22), L * keep, 0.44, fc=CRIMSON if i == 3 else (DARK if i == 0 else tint(0.35 + 0.15 * i)),
                             ec="none"))
    if i: a.add_patch(mp.Rectangle((L * keep, y - 0.22), L * keep / f - L * keep, 0.44, fc="white", ec=GREYF, lw=0.8,
                                   hatch="////", zorder=1))
    a.text(0.12, y, nm, ha="left", va="center", fontsize=8.8, color="white", fontweight="semibold")
    if f:
        a.text(L * keep / f + 0.15, y + 0.07, f"{100 * f:.0f}%", ha="left", va="bottom", fontsize=9.5, color=CRIMSON,
               fontweight="bold")
        a.text(L * keep / f + 0.15, y + 0.02, "lost: " + lost, ha="left", va="top", fontsize=7.6, color=TEXT2)
plt.rcParams["hatch.color"] = GREYF

# (b) sequencing saturation
b = fig.add_subplot(gs[1])
M = 10000                                          # molecules captured from one cell
r = np.linspace(0, 60000, 300); U = M * (1 - np.exp(-r / M))
b.plot(r / 1e3, U / 1e3, color=DARK, lw=2)
b.plot([0, 12], [0, 12], color=LINE, lw=0.9, ls=(0, (3, 2)))
b.text(12.6, 11.4, "if every read found\na new molecule", ha="left", va="center", fontsize=7.6, color=TEXT2,
       linespacing=1.15)
for rr, c in ((5000, TEXT), (40000, CRIMSON)):
    u = M * (1 - np.exp(-rr / M)); sat = 1 - u / rr
    b.scatter([rr / 1e3], [u / 1e3], s=34, color=c, zorder=3)
    b.text(rr / 1e3 + (2.5 if rr < 1e4 else -1), u / 1e3 - (1.2 if rr < 1e4 else 1.5), f"saturation {100 * sat:.0f}%",
           ha="left" if rr < 1e4 else "right", va="top", fontsize=8.4, color=c, fontweight="semibold")
b.text(39, 6.6, "more reads mostly\nre-read old molecules", ha="right", va="top", fontsize=7.6, color=TEXT2,
       linespacing=1.15)
b.set_xlim(0, 60); b.set_ylim(0, 12)
b.set_xlabel("reads per cell (×1,000)", color=TEXT2, fontsize=8.8)
b.set_ylabel("molecules (UMIs) found (×1,000)", color=TEXT2, fontsize=8.8)
for sp in ("left", "bottom"): b.spines[sp].set_color(LINE)
b.tick_params(colors=TEXT2, labelsize=8)
b.set_title("(b) Sequencing saturation", loc="left", fontsize=12, fontweight="semibold", color=TEXT, pad=22)
b.text(0, 1.03, "1 − molecules ÷ reads", transform=b.transAxes, fontsize=8.8, color=TEXT2)

# (c) genes per cell
c = fig.add_subplot(gs[2])
g = np.exp(rng.normal(np.log(2200), 0.45, 4000))
c.hist(g, bins=np.linspace(0, 7000, 50), color=GREYF, ec="white", lw=0.4)
med = np.median(g)
c.axvline(med, color=CRIMSON, lw=1.5)
c.text(med + 150, c.get_ylim()[1] * 0.92, f"median\n{med:,.0f} genes", ha="left", va="top", fontsize=8.4,
       color=CRIMSON, fontweight="semibold", linespacing=1.15)
c.set_xlim(0, 7000); c.set_yticks([])
c.set_xticks([0, 2000, 4000, 6000]); c.set_xticklabels(["0", "2k", "4k", "6k"])
c.set_xlabel("genes detected per cell", color=TEXT2, fontsize=8.8)
for sp in ("bottom",): c.spines[sp].set_color(LINE)
c.spines["left"].set_visible(False)
c.tick_params(colors=TEXT2, labelsize=8)
c.set_title("(c) Median genes per cell", loc="left", fontsize=12, fontweight="semibold", color=TEXT, pad=22)
c.text(0, 1.03, "one number per sample", transform=c.transAxes, fontsize=8.8, color=TEXT2)
fig.text(0.01, 0.02, "Example values from one simulated sample. Compare samples with each other: one sample far below "
         "the rest → look at its dissociation or loading.", fontsize=8.4, color=TEXT2)
pa, pc = a.get_position(), c.get_position()          # rule ends at the right edge of (c), so no blank margin
header(a, "Single-cell RNA-seq · the first QC report", "Five numbers to read before anything else", dy=40,
       x1=(pc.x1 - pa.x0) / pa.width)
fig.savefig(OUT + "qc_metrics.png"); print("ok k33 qc metrics", round(med))
