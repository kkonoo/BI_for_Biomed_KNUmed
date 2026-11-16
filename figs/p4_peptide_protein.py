from gallery import *
import os
rng = np.random.default_rng(11)
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "assets", "img", "w6") + "/"

# Peptide → protein (simulated): five peptides of one protein ionise differently, so they sit at different
# heights, but they move together between groups. A summary (Tukey median polish, as in MSstats) removes the
# peptide effect and gives one value per sample. The weakest peptide drops out below the detection limit.
fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.2), gridspec_kw={"wspace": 0.32, "width_ratios": [1.35, 1]})
DARK = "#5E6670"; GREYF = "#B8BEC7"

def style(a, title, sub):
    for sp in ("left", "bottom"): a.spines[sp].set_color(LINE)
    a.tick_params(colors=TEXT2, labelsize=8.3)
    a.set_title(title, loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
    a.text(0, 1.03, sub, transform=a.transAxes, fontsize=9, color=TEXT2)

base = np.array([26.2, 24.8, 23.6, 22.5, 21.1])
grp = np.array([0, 0, 0, 1, 1, 1]); eff = 1.0
Y = base[:, None] + eff * grp[None, :] + rng.normal(0, 0.15, (5, 6)) + rng.normal(0, 0.12, 6)[None, :]
Y[4, [0, 2]] = np.nan                       # the weakest peptide is not seen in two control samples
x = np.arange(6)

a = axes[0]
a.axvspan(2.5, 5.5, color=FAINT, zorder=0)
for k in range(5):
    col = [DARK, DARK, DARK, DARK, GREYF][k]
    a.plot(x, Y[k], color=col, lw=1.1, zorder=2)
    a.scatter(x, Y[k], s=26, color=col, zorder=3, lw=0)
    a.text(5.25, Y[k, 5], f"peptide {k + 1}", fontsize=8.4, color=TEXT2, va="center", ha="left")
a.scatter([0, 2], [20.2, 20.2], s=26, facecolor="white", edgecolor=CRIMSON, lw=1.1, zorder=3)
a.axhline(20.6, color=TEXT2, lw=0.8, ls=(0, (3, 2)))
a.text(6.55, 20.7, "detection limit", fontsize=8.2, color=TEXT2, va="bottom", ha="right")
a.text(0.1, 19.75, "missing: too faint to see", fontsize=8.2, color=CRIMSON, va="top", ha="left")
a.set_xticks(x); a.set_xticklabels(["C1", "C2", "C3", "T1", "T2", "T3"])
a.set_xlim(-0.4, 6.6); a.set_ylim(19.2, 28.2)
a.set_ylabel("log₂ intensity", color=TEXT2, fontsize=9)
a.text(-0.3, 27.45, "height = how well each peptide ionises,\nnot how much protein there is", fontsize=8.3,
       color=CRIMSON, ha="left", va="center")
style(a, "(a) Five peptides of one protein", "C = control, T = treated; they rise together")
a.spines["bottom"].set_bounds(-0.4, 5.4)

# (b) Tukey median polish
def median_polish(Y, it=10):
    R = Y.copy(); mu = 0.0; r = np.zeros(Y.shape[0]); c = np.zeros(Y.shape[1])
    for _ in range(it):
        d = np.nanmedian(R, 1); r += d; R -= d[:, None]
        d = np.nanmedian(c); mu += d; c -= d
        d = np.nanmedian(R, 0); c += d; R -= d[None, :]
        d = np.nanmedian(r); mu += d; r -= d
    return mu, r, c
mu, r, cc = median_polish(Y)
prot = mu + cc
b = axes[1]
b.axvspan(2.5, 5.5, color=FAINT, zorder=0)
b.scatter(x, prot, s=46, color=CRIMSON, zorder=3, lw=0)
m0, m1 = prot[:3].mean(), prot[3:].mean()
b.plot([-0.2, 2.2], [m0, m0], color=TEXT2, lw=1); b.plot([2.8, 5.2], [m1, m1], color=TEXT2, lw=1)
b.annotate("", (3.0, m1), (3.0, m0), arrowprops=dict(arrowstyle="-|>", color=TEXT, lw=1, mutation_scale=10,
                                                     shrinkA=0, shrinkB=0))
b.text(2.6, (m0 + m1) / 2, f"log₂FC\n= {m1 - m0:.2f}", fontsize=9, color=TEXT, ha="right", va="center")
b.set_xticks(x); b.set_xticklabels(["C1", "C2", "C3", "T1", "T2", "T3"])
b.set_xlim(-0.4, 5.4); b.set_ylim(m0 - 1.2, m1 + 1.2)
b.set_ylabel("protein estimate (log₂)", color=TEXT2, fontsize=9)
b.text(0.02, 0.04, "peptide effect removed;\nsamples with a missing peptide\nare still estimated from the rest",
       transform=b.transAxes, fontsize=8.3, color=TEXT2, va="bottom")
style(b, "(b) One value per sample", "Tukey median polish (MSstats)")

header(a, "Proteomics · peptide to protein", "A protein value is a summary of its peptides", dy=40, x1=2.2)
fig.savefig(OUT + "peptide_to_protein.png"); print("ok p4", round(m1 - m0, 2))
