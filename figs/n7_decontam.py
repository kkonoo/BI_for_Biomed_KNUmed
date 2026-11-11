from gallery import *
rng = np.random.default_rng(7)

# How decontam-type tools spot kit DNA (simulated). Each extraction adds the same small amount of kit DNA,
# so (a) a kit taxon's share falls as 1 / total DNA, while a real taxon's share does not depend on it;
# (b) a kit taxon makes up a large share of the blank extractions, a real one almost none.
n, nb = 60, 12
sample_dna = np.exp(rng.uniform(np.log(0.05), np.log(50), n))           # ng of bacterial DNA from the sample
kit_dna = 0.04 * np.exp(rng.normal(0, 0.3, n + nb))                     # ng from the kit, every extraction
real_share = 0.08 * np.exp(rng.normal(0, 0.35, n))                      # one genuine taxon, share of sample DNA
kit_frac = 0.05                                                          # one kit genus, share of kit DNA
tot = sample_dna + kit_dna[:n]
f_real = real_share * sample_dna / tot
f_kit = kit_frac * kit_dna[:n] / tot * np.exp(rng.normal(0, 0.25, n))

fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.3), gridspec_kw={"wspace": 0.3, "width_ratios": [1.2, 1]})
DARK = "#5E6670"; GREYF = "#B8BEC7"
a = axes[0]
a.scatter(tot, f_real, s=16, color=DARK, lw=0, label="a real gut taxon")
a.scatter(tot, f_kit, s=16, color=CRIMSON, lw=0, label="a kit genus")
xx = np.logspace(np.log10(tot.min()), np.log10(tot.max()), 50)
a.plot(xx, kit_frac * 0.04 / xx, color=CRIMSON, lw=0.9, ls=(0, (3, 2)))
a.plot(xx, np.full_like(xx, np.median(f_real)), color=DARK, lw=0.9, ls=(0, (3, 2)))
a.set_xscale("log"); a.set_yscale("log")
a.set_xlabel("total DNA in the extract (ng)", color=TEXT2, fontsize=9)
a.set_ylabel("share of reads", color=TEXT2, fontsize=9)
a.set_yticks([1e-4, 1e-3, 1e-2, 1e-1]); a.set_yticklabels(["0.01%", "0.1%", "1%", "10%"])
a.set_ylim(top=0.6)
a.text(xx[-1], 0.3, "flat: comes from the sample", ha="right", va="center", fontsize=8.5,
       color=DARK, fontweight="semibold")
a.text(1.3, kit_frac * 0.04 / 1.3 * 1.8, "falls as 1 / DNA:\nsame amount added every time", ha="left", va="bottom",
       fontsize=8.5, color=CRIMSON, fontweight="semibold", linespacing=1.15)
for sp in ("top", "right"): a.spines[sp].set_visible(False)
for sp in ("left", "bottom"): a.spines[sp].set_color(LINE)
a.tick_params(colors=TEXT2, labelsize=8.3); a.minorticks_off()
a.set_title("(a) Frequency vs DNA amount", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
a.text(0, 1.03, "60 samples with very different biomass", transform=a.transAxes, fontsize=9, color=TEXT2)

b = axes[1]
# share of reads: median over the 60 samples vs over the 12 blanks (blanks contain only kit DNA,
# plus a trace of the real taxon from cross-contamination)
blank_real_share = 3e-6 * np.exp(rng.normal(0, 0.5, nb))
vals = [("a real gut taxon", np.median(f_real), np.median(blank_real_share), DARK),
        ("a kit genus", np.median(f_kit), kit_frac, CRIMSON)]
fmt = lambda v: f"{100 * v:.2g}%" if v >= 1e-4 else f"{100 * v:.1g}%"
for k, (nm, vs, vb, col) in enumerate(vals):
    b.bar(k - 0.18, vs, width=0.32, color=GREYF)
    b.bar(k + 0.18, vb, width=0.32, color=col)
    b.text(k - 0.18, vs * 1.25, fmt(vs), ha="center", va="bottom", fontsize=8.6, color=TEXT)
    b.text(k + 0.18, vb * 1.25, fmt(vb), ha="center", va="bottom", fontsize=8.6, color=TEXT, fontweight="semibold")
b.set_yscale("log"); b.set_ylim(1e-6, 1)
b.set_yticks([1e-6, 1e-4, 1e-2, 1]); b.set_yticklabels(["0.0001%", "0.01%", "1%", "100%"])
b.set_xticks([0, 1]); b.set_xticklabels([v[0] for v in vals], fontsize=9, color=TEXT)
for sp in ("top", "right"): b.spines[sp].set_visible(False)
for sp in ("left", "bottom"): b.spines[sp].set_color(LINE)
b.tick_params(colors=TEXT2, labelsize=8.3); b.tick_params(axis="x", length=0); b.minorticks_off()
b.set_ylabel("share of reads (median)", color=TEXT2, fontsize=9)
b.set_title("(b) Samples vs blanks", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
b.text(0, 1.03, f"grey = {n} samples, coloured = {nb} blank extractions", transform=b.transAxes, fontsize=9,
       color=TEXT2)
header(a, "Metagenomics · decontamination", "Kit DNA leaves two fingerprints: it tracks low DNA, and it is in the blanks",
       dy=40, x1=2.25)
fig.savefig(OUT + "decontam.png")
print("ok n7", [(v[0], v[1], v[2]) for v in vals])
