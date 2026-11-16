from gallery import *
import os
rng = np.random.default_rng(6)
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "assets", "img", "w6") + "/"

# Proteomics QC, one panel per step (illustrative numbers):
# (a) identification — proteins per run, one run far below the rest;
# (b) contaminants — keratin share of intensity, one sample high;
# (c) quantification — CV of a pooled QC sample (technical) vs. replicates within a group (technical + biological).
fig, axes = plt.subplots(1, 3, figsize=(10.4, 4.1), gridspec_kw={"wspace": 0.34, "width_ratios": [1.1, 1.1, 1]})
DARK = "#5E6670"; GREYF = "#B8BEC7"

def style(a, title, sub):
    for sp in ("left", "bottom"): a.spines[sp].set_color(LINE)
    a.tick_params(colors=TEXT2, labelsize=8.3)
    a.set_title(title, loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
    a.text(0, 1.03, sub, transform=a.transAxes, fontsize=9, color=TEXT2)

# (a) proteins identified per run
a = axes[0]
ids = np.round(rng.normal(5200, 120, 12)).astype(int); ids[6] = 2650
x = np.arange(1, 13)
a.plot(x, ids, color=LINE, lw=1.0, zorder=1)
a.scatter(x, ids, s=34, color=[CRIMSON if i == 6 else DARK for i in range(12)], zorder=3, lw=0)
med = np.median(ids)
a.axhline(med, color=TEXT2, lw=0.8, ls=(0, (3, 2)), zorder=0)
a.text(12.4, 5650, "dashed = median of all runs", fontsize=8.2, color=TEXT2, ha="right", va="bottom")
a.annotate("run 7: half the proteins\n→ find the cause (LC column,\ninjection) before using it",
           (7, ids[6]), (1.0, 1300), fontsize=8.4, color=CRIMSON, ha="left", va="center",
           arrowprops=dict(arrowstyle="-", color=CRIMSON, lw=0.7))
a.set_xlim(0.3, 12.7)
a.set_xticks([1, 4, 7, 10, 12]); a.set_xlabel("run (injection order)", color=TEXT2, fontsize=9)
a.set_ylabel("proteins identified (1% FDR)", color=TEXT2, fontsize=9); a.set_ylim(0, 6300)
style(a, "(a) Identification", "proteins per run")

# (b) keratin share of total intensity
b = axes[1]
ker = np.clip(rng.normal(0.35, 0.08, 12), 0.15, None); ker[3] = 4.8
cols = [GREYF] * 12; cols[3] = CRIMSON
b.bar(np.arange(1, 13), ker, color=cols, width=0.7)
b.annotate("sample 4: skin and hair keratin\n→ a handling problem,\nnot biology", (4, ker[3]), (5.2, 3.9),
           fontsize=8.4, color=CRIMSON, ha="left", va="center",
           arrowprops=dict(arrowstyle="-", color=CRIMSON, lw=0.7))
b.set_xticks([1, 4, 7, 10, 12]); b.set_xlabel("sample", color=TEXT2, fontsize=9)
b.set_ylabel("keratin share of intensity (%)", color=TEXT2, fontsize=9); b.set_ylim(0, 5.6)
style(b, "(b) Contaminants", "keratin, the classic contaminant")

# (c) CV distributions
c = axes[2]
cv_qc = np.exp(rng.normal(np.log(8), 0.45, 3000))
cv_bio = np.exp(rng.normal(np.log(24), 0.45, 3000))
bins = np.linspace(0, 70, 57)
c.hist(cv_qc, bins=bins, color=DARK, alpha=0.9, lw=0)
c.hist(cv_bio, bins=bins, color=GREYF, alpha=0.9, lw=0)
c.annotate(f"pooled QC sample,\ninjected 8 times\nmedian CV {np.median(cv_qc):.0f}%", (9.5, 430), (20, 600),
           fontsize=8.3, color=DARK, ha="left", va="center",
           arrowprops=dict(arrowstyle="-", color=DARK, lw=0.7))
c.annotate(f"3 replicates in one group\nmedian CV {np.median(cv_bio):.0f}%", (30, 150), (40, 250),
           fontsize=8.3, color=TEXT2, ha="left", va="center",
           arrowprops=dict(arrowstyle="-", color=TEXT2, lw=0.7))
c.text(69, 420, "technical noise\nmust be well below\nthe biological spread", fontsize=8.4, color=CRIMSON,
       ha="right", va="bottom")
c.set_xlim(0, 70); c.set_ylim(0, 700); c.set_yticks([]); c.spines["left"].set_visible(False)
c.set_xlabel("CV per protein (%)", color=TEXT2, fontsize=9)
style(c, "(c) Quantification", "CV = SD ÷ mean, per protein")
c.spines["left"].set_visible(False)

header(a, "Proteomics · QC", "Check each step against the other runs, not against a fixed number",
       dy=40, x1=3.25)
fig.savefig(OUT + "prot_qc.png"); print("ok p2 qc", med, np.median(cv_qc), np.median(cv_bio))
