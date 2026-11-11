from gallery import *
from scipy import stats
rng = np.random.default_rng(11)

# TMT batch structure: 90 samples (45 patients, 45 controls) in 18-plex sets.
# (a) confounded layout; (b) randomised layout with a shared reference channel; (c) false-positive rate
# for 2,000 proteins with batch effects and no true difference (simulated).
fig = plt.figure(figsize=(9.8, 4.4))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.72], wspace=0.22, left=0.03, right=0.98, top=0.78, bottom=0.12)
DARK = "#5E6670"; GREYF = "#B8BEC7"
PAT, CON, REF = CRIMSON, GREYF, DARK

# layouts: lists of batches, each a list of channel labels
conf = []
labels = ["P"] * 45 + ["C"] * 45
for k in range(5): conf.append(labels[k * 18:(k + 1) * 18])
bal = []
for k, npat in enumerate((8, 8, 8, 7, 7, 7)):                 # stratified: every set gets 7–8 of each group
    bal.append(["R"] + list(rng.permutation(["P"] * npat + ["C"] * (15 - npat))) + ["-", "-"])

def grid(ax, batches, title, sub):
    ax.set_xlim(-2.2, 18.2); ax.set_ylim(-0.6, 6.4); ax.axis("off")
    for r, bt in enumerate(batches):
        y = 5.4 - r
        ax.text(-0.4, y + 0.4, f"set {r + 1}", ha="right", va="center", fontsize=8, color=TEXT2)
        for j, c in enumerate(bt):
            if c == "-": continue
            fc = {"P": PAT, "C": CON, "R": REF}[c]
            ax.add_patch(mp.Rectangle((j + 0.06, y + 0.08), 0.88, 0.72, fc=fc, ec="none"))
    ax.set_title(title, loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
    ax.text(0, 1.03, sub, transform=ax.transAxes, fontsize=9, color=TEXT2)

a = fig.add_subplot(gs[0]); grid(a, conf, "(a) Confounded", "patients in sets 1–3, controls in 3–5")
b = fig.add_subplot(gs[1]); grid(b, bal, "(b) Randomised + reference", "each set: a mix of both, plus one shared channel")
b.text(0.5, -0.45, "dark = the same reference sample in every set", ha="left", va="center", fontsize=8.3, color=TEXT2)
a.text(0, -0.45, "red = patient   grey = control", ha="left", va="center", fontsize=8.3, color=TEXT2)

# (c) simulation: 2,000 proteins, batch effects sd 0.5, sample sd 0.3, no group difference
NP = 2000
def fpr(batches, use_ref):
    hits = 0
    for _ in range(NP):
        be = rng.normal(0, 0.5, len(batches))
        g1, g2 = [], []
        for k, bt in enumerate(batches):
            ref = be[k] + rng.normal(0, 0.1)
            for c in bt:
                if c in "PC":
                    v = be[k] + rng.normal(0, 0.3) - (ref if use_ref else 0)
                    (g1 if c == "P" else g2).append(v)
        hits += stats.ttest_ind(g1, g2).pvalue < 0.05
    return 100 * hits / NP
f_conf, f_bal = fpr(conf, False), fpr(bal, True)
c = fig.add_subplot(gs[2])
c.bar([0, 1], [f_conf, f_bal], color=[CRIMSON, DARK], width=0.55)
for x, v in zip([0, 1], [f_conf, f_bal]):
    c.text(x, v + 2, f"{v:.0f}%", ha="center", va="bottom", fontsize=10.5, fontweight="bold", color=TEXT)
c.axhline(5, color=TEXT2, lw=0.8, ls=(0, (3, 2)))
c.set_xticks([0, 1]); c.set_xticklabels(["(a)", "(b)"], fontsize=9.5, color=TEXT)
c.tick_params(axis="x", length=0); c.set_yticks([]); c.set_xlim(-0.5, 1.7); c.set_ylim(0, 105)
for sp in ("left", "top", "right"): c.spines[sp].set_visible(False)
c.spines["bottom"].set_color(LINE)
c.set_title("(c) 'Significant' proteins", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
c.text(0, 1.03, "no true difference; dashed = 5% by chance", transform=c.transAxes, fontsize=9, color=TEXT2)
header(a, "Proteomics · TMT batches", "90 samples do not fit in one TMT set — plan the sets", dy=40, x1=3.1)
fig.savefig(OUT + "tmt_batches.png"); print("ok m11", round(f_conf), round(f_bal))
