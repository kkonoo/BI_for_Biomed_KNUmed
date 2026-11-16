from gallery import *
from scipy import stats
import os
rng = np.random.default_rng(5)
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "assets", "img", "w6") + "/"

# Enrichment background in proteomics (simulated): MS detects abundant, housekeeping proteins far more often than
# receptors or transcription factors. A "hit list" of 300 proteins drawn at random from the detected proteins —
# no biology at all — looks enriched for ribosomes and glycolysis when tested against the whole genome,
# and shows nothing when tested against the proteins that could have been detected.
cats = [("Ribosome", 150, 0.95), ("Glycolysis / TCA", 110, 0.9), ("Proteasome", 45, 0.85),
        ("Spliceosome", 130, 0.75), ("Transcription factors", 1600, 0.10), ("GPCRs", 800, 0.03),
        ("Olfactory receptors", 400, 0.0)]
G = 20000
sizes = np.array([c[1] for c in cats]); rates = np.array([c[2] for c in cats])
det_in = np.round(sizes * rates).astype(int)
D = 5000                                         # proteins detected in the experiment
rest_det = D - det_in.sum()
# random hit list of 300 drawn from the detected proteins
labels = np.concatenate([np.full(n, i) for i, n in enumerate(det_in)] + [np.full(rest_det, -1)])
hits = rng.choice(labels, 300, replace=False)
k = np.array([(hits == i).sum() for i in range(len(cats))])
p_gen = stats.hypergeom.sf(k - 1, G, sizes, 300)
p_det = stats.hypergeom.sf(k - 1, D, det_in, 300)

fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.1), gridspec_kw={"wspace": 0.75, "width_ratios": [1, 1.05]})
DARK = "#5E6670"; GREYF = "#B8BEC7"
y = np.arange(len(cats))[::-1]

def style(a, title, sub):
    for sp in ("left", "bottom"): a.spines[sp].set_color(LINE)
    a.tick_params(colors=TEXT2, labelsize=8.3)
    a.set_title(title, loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
    a.text(0, 1.03, sub, transform=a.transAxes, fontsize=9, color=TEXT2)

a = axes[0]
a.barh(y, 100 * rates, color=[DARK if r > 0.5 else GREYF for r in rates], height=0.62)
for yy, r in zip(y, rates):
    a.text(100 * r + 2, yy, f"{100 * r:.0f}%", va="center", fontsize=8.6, color=TEXT)
a.set_yticks(y); a.set_yticklabels([c[0] for c in cats], fontsize=9, color=TEXT)
a.tick_params(axis="y", length=0); a.set_xlim(0, 112)
a.set_xlabel("% of the gene set detected by MS", color=TEXT2, fontsize=9)
style(a, "(a) What MS can see", "abundant housekeeping proteins first")

b = axes[1]
lg, ld = -np.log10(p_gen), -np.log10(p_det)
b.barh(y + 0.17, lg, height=0.32, color=CRIMSON, label="background = whole genome (20,000)")
b.barh(y - 0.17, ld, height=0.32, color=GREYF, label="background = detected proteins (5,000)")
b.axvline(-np.log10(0.05), color=TEXT2, lw=0.8, ls=(0, (3, 2)))
b.text(-np.log10(0.05) - 0.08, y[-1] - 0.75, "p = 0.05", fontsize=8.2, color=TEXT2, va="center", ha="right")
b.set_yticks(y); b.set_yticklabels([c[0] for c in cats], fontsize=9, color=TEXT)
b.tick_params(axis="y", length=0)
b.set_xlabel("enrichment, −log₁₀ p", color=TEXT2, fontsize=9)
b.legend(frameon=True, framealpha=1, edgecolor="none", fontsize=8.4, loc="center right", bbox_to_anchor=(1.02, 0.36), labelcolor=TEXT, handlelength=1.2)
b.set_ylim(-1.1, len(cats) - 0.4)
style(b, "(b) A random list of 300 proteins", "no biology in it at all")

header(a, "Proteomics · enrichment", "Test against what you could have detected, not the whole genome",
       dy=40, x1=2.95)
fig.savefig(OUT + "prot_background.png"); print("ok p5", k, np.round(lg, 1), np.round(ld, 1))
