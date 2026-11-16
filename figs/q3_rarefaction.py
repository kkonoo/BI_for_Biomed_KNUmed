from gallery import *
import os
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "assets", "img", "w7") + "/"
rng = np.random.default_rng(5)

# Rarefaction (simulated communities): (a) expected observed richness when n reads are drawn,
# E[S_n] = Σ 1 − (1 − p_i)^n, for samples of different depth and richness; the cut at 8,000 reads.
# (b) what rarefying to that depth does to each library: reads thrown away, and a sample dropped.
DARK = "#5E6670"; GREYF = "#B8BEC7"
def community(S, sd):
    p = np.exp(rng.normal(0, sd, S)); return p / p.sum()
# (name, number of taxa, spread of abundances, library depth); a wide spread = many rare taxa
samples = [("A", 250, 1.6, 52000), ("B", 400, 2.0, 38000), ("C", 150, 1.0, 24000), ("D", 900, 2.6, 61000),
           ("E", 200, 1.3, 12000), ("F", 700, 2.4, 3000)]
CUT = 8000

fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.3), gridspec_kw={"wspace": 0.32, "width_ratios": [1.35, 1]})
def tidy(a):
    for sp in ("top", "right"): a.spines[sp].set_visible(False)
    for sp in ("left", "bottom"): a.spines[sp].set_color(LINE)
    a.tick_params(colors=TEXT2, labelsize=8.3); a.minorticks_off()

a = axes[0]
a.axvspan(CUT, 66000, color=FAINT, zorder=0)
for name, S, sd, depth in samples:
    p = community(S, sd)
    n = np.linspace(1, depth, 200)
    rich = np.array([(1 - (1 - p) ** k).sum() for k in n])
    col = CRIMSON if name == "D" else (DARK if name == "C" else GREYF)
    lw = 1.8 if name in "CD" else 1.3
    a.plot(n, rich, color=col, lw=lw, zorder=3)
    a.scatter([depth], [rich[-1]], s=16, color=col, zorder=4, lw=0)
    a.text(depth + 900, rich[-1], name, ha="left", va="center", fontsize=9, color=col, fontweight="semibold")
    if name == "D": rD = rich
    if name == "C": rC = rich
a.axvline(CUT, color=TEXT, lw=1, ls=(0, (3, 2)), zorder=2)
a.text(CUT + 700, 8, "rarefy every sample\nto 8,000 reads", ha="left", va="bottom", fontsize=8.6, color=TEXT,
       linespacing=1.15)
a.text(50000, 520, "D: still rising →\nmore taxa left to find", ha="center", va="top", fontsize=8.6,
       color=CRIMSON, fontweight="semibold", linespacing=1.15)
a.text(30000, rC[-1] - 14, "C: levels off →\ndepth was enough", ha="center", va="top", fontsize=8.6,
       color=DARK, fontweight="semibold", linespacing=1.15)
a.text(65000, 12, "shaded: reads beyond the cut,\nnot used after rarefying", ha="right", va="bottom", linespacing=1.15, fontsize=8.4, color=TEXT2)
a.set_xlim(0, 66000); a.set_ylim(0, 700)
a.set_xticks([0, 20000, 40000, 60000]); a.set_xticklabels(["0", "20k", "40k", "60k"])
a.set_xlabel("reads drawn from the library", color=TEXT2, fontsize=9)
a.set_ylabel("taxa observed", color=TEXT2, fontsize=9)
tidy(a)
a.set_title("(a) Rarefaction curves", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
a.text(0, 1.03, "six samples; dot = the library's full depth", transform=a.transAxes, fontsize=9, color=TEXT2)

b = axes[1]
yy = np.arange(len(samples))[::-1]
for yv, (name, S, sd, depth) in zip(yy, samples):
    if depth >= CUT:
        b.barh(yv, CUT, color=CRIMSON, height=0.6)
        b.barh(yv, depth - CUT, left=CUT, color=GREYF, height=0.6)
        b.text(depth + 1200, yv, f"{100 * (depth - CUT) / depth:.0f}% thrown away", ha="left", va="center",
               fontsize=8.4, color=TEXT2)
    else:
        b.barh(yv, depth, color="white", edgecolor=GREYF, height=0.6, hatch="////", lw=0.8)
        b.text(depth + 1200, yv, "below the cut → sample dropped", ha="left", va="center", fontsize=8.4,
               color=CRIMSON, fontweight="semibold")
b.axvline(CUT, color=TEXT, lw=1, ls=(0, (3, 2)))
b.set_yticks(yy); b.set_yticklabels([s[0] for s in samples], fontsize=9, color=TEXT)
b.set_xlim(0, 90000); b.set_xticks([0, 20000, 40000, 60000]); b.set_xticklabels(["0", "20k", "40k", "60k"])
b.set_xlabel("reads in the library", color=TEXT2, fontsize=9)
tidy(b); b.tick_params(axis="y", length=0)
b.set_title("(b) What rarefying costs", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
b.text(0, 1.03, "red = reads kept (8,000 each), grey = discarded", transform=b.transAxes, fontsize=9,
       color=TEXT2)
header(a, "Metagenomics · sequencing depth", "Richness grows with depth — rarefying evens it out by throwing reads away",
       dy=40, x1=2.1)
fig.savefig(OUT + "rarefaction.png"); print("ok q3 rarefaction")
