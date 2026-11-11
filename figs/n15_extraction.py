from gallery import *

# Extraction bias: a mock community of eight species in equal amounts, extracted two ways. Gram-positive cell walls
# are tough; without bead-beating they break open poorly. Lysis efficiencies are illustrative.
species = [("Escherichia coli", "−", 0.95, 0.95), ("Pseudomonas aeruginosa", "−", 0.95, 0.95),
           ("Salmonella enterica", "−", 0.95, 0.95),
           ("Bacillus subtilis", "+", 0.80, 0.12), ("Enterococcus faecalis", "+", 0.85, 0.30),
           ("Staphylococcus aureus", "+", 0.85, 0.20), ("Listeria monocytogenes", "+", 0.90, 0.35),
           ("Limosilactobacillus fermentum", "+", 0.85, 0.25)]
bead = np.array([s[2] for s in species]); gentle = np.array([s[3] for s in species])
bead /= bead.sum(); gentle /= gentle.sum()
fig, ax = plt.subplots(figsize=(9.8, 4.4))
fig.subplots_adjust(left=0.32, right=0.8)
DARK = "#5E6670"; GREYF = "#B8BEC7"
Y = np.arange(len(species))[::-1]
ax.axvline(1 / 8, color=TEXT, lw=0.8, ls=(0, (3, 2)))
ax.text(1 / 8 + 0.003, len(species) - 0.25, "true share 12.5%", ha="left", va="bottom", fontsize=8.6, color=TEXT)
for y, b_, g_ in zip(Y, bead, gentle):
    ax.plot([b_, g_], [y, y], color=LINE, lw=1.2, zorder=1)
ax.scatter(bead, Y, s=48, color=DARK, zorder=3, label="with bead-beating")
ax.scatter(gentle, Y, s=48, color=CRIMSON, zorder=3, label="gentle lysis only")
ax.set_yticks(Y); ax.set_yticklabels([s[0] for s in species], fontsize=9, color=TEXT)
ax.tick_params(axis="y", length=0); ax.tick_params(colors=TEXT2, labelsize=8.3)
for lo, hi in ((Y[2], Y[0]), (Y[-1], Y[3])):
    ax.plot([1.01] * 2, [lo - 0.3, hi + 0.3], color=LINE, lw=1.2, transform=ax.get_yaxis_transform(), clip_on=False)
ax.text(1.03, (Y[0] + Y[2]) / 2, "Gram-negative:\nthin wall, easy to open", transform=ax.get_yaxis_transform(),
        ha="left", va="center", fontsize=8.6, color=TEXT, linespacing=1.15)
ax.text(1.03, (Y[3] + Y[-1]) / 2, "Gram-positive:\nthick wall; without\nbead-beating, under-counted",
        transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=8.6, color=TEXT, linespacing=1.15)
ax.set_xlim(0, 0.26); ax.set_ylim(-0.6, len(species) + 0.2)
ax.set_xticks([0, 0.05, 0.1, 0.15, 0.2, 0.25]); ax.set_xticklabels(["0", "5%", "10%", "15%", "20%", "25%"])
ax.set_xlabel("share of reads", color=TEXT2, fontsize=9)
for sp in ("top", "right", "left"): ax.spines[sp].set_visible(False)
ax.spines["bottom"].set_color(LINE)
ax.legend(frameon=False, fontsize=8.6, loc="lower right", labelcolor=TEXT, handletextpad=0.2)
XL = -0.6                                           # align titles with the species names, left of the axes
ax.text(XL, 1.03, "a mock community; lysis efficiencies illustrative", transform=ax.transAxes, fontsize=9, color=TEXT2)
ax.set_title("Eight species, equal amounts, two protocols", loc="left", x=XL, fontsize=12.3, fontweight="semibold",
             color=TEXT, pad=22)
header(ax, "Metagenomics · DNA extraction", "The extraction protocol decides which bacteria you see", x=XL, dy=40,
       x1=1.45)
fig.savefig(OUT + "extraction_bias.png"); print("ok n15", np.round(bead, 3), np.round(gentle, 3))
