from gallery import *

# Poly(A) selection vs rRNA depletion: where the reads go, how each works, and when to use it
fig = plt.figure(figsize=(9.8, 5.6))
gs = fig.add_gridspec(1, 2, width_ratios=[0.8, 1.45], wspace=0.12, left=0.06, right=0.995, top=0.8, bottom=0.2)
DARK = "#5E6670"; GREYF = "#B8BEC7"

# (a) where the reads go -----------------------------------------------------------------------------
a = fig.add_subplot(gs[0])
cats = [("rRNA", DARK), ("mRNA", CRIMSON), ("lncRNA / other non-coding", tint(0.45)), ("pre-mRNA (introns)", tint(0.18))]
vals = {"total RNA\n(no selection)": [0.85, 0.03, 0.10, 0.02],
        "poly(A)\nselection": [0.02, 0.86, 0.09, 0.03],
        "rRNA\ndepletion": [0.05, 0.55, 0.29, 0.11]}
X = np.arange(3); bottom = np.zeros(3)
for i, (lab, col) in enumerate(cats):
    v = np.array([vals[k][i] for k in vals])
    a.bar(X, v, bottom=bottom, color=col, width=0.62, edgecolor="white", lw=1.2, label=lab)
    bottom += v
a.set_xticks(X); a.set_xticklabels(list(vals), fontsize=8.8, color=TEXT)
a.set_ylim(0, 1); a.set_yticks([0, 0.5, 1]); a.set_yticklabels(["0", "50%", "100%"])
a.set_ylabel("share of sequenced reads", color=TEXT2, fontsize=9)
for sp in ("left", "bottom"): a.spines[sp].set_color(LINE)
a.tick_params(colors=TEXT2, labelsize=8.5)
a.text(0, 0.45, "85%\nrRNA", ha="center", va="center", fontsize=9, color="white", fontweight="bold", linespacing=1.2)
a.legend(frameon=False, fontsize=8.4, loc="upper center", bbox_to_anchor=(0.45, -0.17), ncol=2,
         handlelength=1, handletextpad=0.4, columnspacing=0.8, labelcolor=TEXT)
a.set_title("(a) Where the reads go", loc="left", fontsize=12.5, fontweight="semibold", color=TEXT, pad=22)
a.text(0, 1.03, "without selection, rRNA eats the budget", transform=a.transAxes, ha="left", va="bottom",
       fontsize=9.3, color=TEXT2)

# (b) how each works and when to use it -----------------------------------------------------------------
b = fig.add_subplot(gs[1]); b.set_xlim(0, 10); b.set_ylim(0, 10); b.axis("off")
b.set_title("(b) Two ways to get rid of rRNA", loc="left", fontsize=12.5, fontweight="semibold", color=TEXT, pad=22)
b.text(0, 1.03, "decided at the bench — cannot be undone in analysis", transform=b.transAxes, ha="left",
       va="bottom", fontsize=9.3, color=TEXT2)
def card(x0, title, pros, cons, use):
    b.add_patch(mp.FancyBboxPatch((x0, 0.0), 4.75, 10, boxstyle="round,pad=0,rounding_size=0.25",
                                  fc="#F8F9FA", ec=LINE, lw=0.9))
    b.text(x0 + 0.3, 9.45, title, ha="left", va="center", fontsize=11.5, fontweight="bold", color=TEXT)
    y = 5.7
    for t in pros:
        b.text(x0 + 0.3, y, "✓", ha="left", va="center", fontsize=9.5, color=CRIMSON, fontweight="bold")
        b.text(x0 + 0.75, y, t, ha="left", va="center", fontsize=9.2, color=TEXT); y -= 0.78
    for t in cons:
        b.text(x0 + 0.3, y, "✕", ha="left", va="center", fontsize=9.5, color=DARK, fontweight="bold")
        b.text(x0 + 0.75, y, t, ha="left", va="center", fontsize=9.2, color=TEXT2); y -= 0.78
    b.text(x0 + 0.3, 0.75, use, ha="left", va="center", fontsize=9, color=CRIMSON, fontweight="semibold")
# poly(A): oligo-dT bead catches tails
card(0.0, "poly(A) selection",
     ["clean, cheap per usable read", "the default for DE studies"],
     ["needs intact RNA (fails on FFPE)", "misses RNA without a tail", "3′ bias if RNA is degraded"],
     "→ standard mRNA DE, good RNA")
b.add_patch(mp.Circle((1.0, 7.6), 0.55, fc=FAINT, ec=LINE, lw=1))
b.text(1.0, 7.6, "bead", ha="center", va="center", fontsize=7.5, color=TEXT2)
b.plot([1.55, 2.2], [7.6, 7.6], color=TEXT2, lw=1.4)
b.text(1.88, 7.25, "TTTT", ha="center", va="center", fontsize=7.5, color=TEXT2, family="DejaVu Sans Mono")
b.plot([2.25, 4.3], [7.95, 7.95], color=CRIMSON, lw=2.4)
b.text(2.35, 7.6, "AAAA", ha="left", va="center", fontsize=7.5, color=CRIMSON, family="DejaVu Sans Mono")
b.text(3.3, 8.3, "mRNA", ha="center", va="center", fontsize=8.5, color=CRIMSON)
b.plot([2.4, 4.3], [6.7, 6.7], color=DARK, lw=2.4)
b.text(3.35, 6.4, "rRNA (no tail) washed away", ha="center", va="center", fontsize=7.8, color=TEXT2)
# depletion: probes pull rRNA out
card(5.25, "rRNA depletion",
     ["works on degraded / FFPE RNA", "keeps non-coding and nascent RNA"],
     ["more reads for the same mRNA depth", "costs more; intronic background"],
     "→ FFPE, lncRNA, circRNA, nascent RNA")
b.plot([5.75, 7.55], [8.0, 8.0], color=DARK, lw=2.4)
for xx in np.arange(5.85, 7.5, 0.3):
    b.plot([xx, xx], [7.72, 7.92], color=LINE, lw=1)
b.plot([5.75, 7.55], [7.62, 7.62], color=TEXT2, lw=1.4)
b.text(6.65, 8.3, "rRNA", ha="center", va="center", fontsize=8.5, color=DARK)
b.text(6.65, 7.25, "probe binds → removed", ha="center", va="center", fontsize=7.8, color=TEXT2)
b.plot([7.95, 9.75], [8.0, 8.0], color=CRIMSON, lw=2.4)
b.text(8.85, 8.3, "mRNA", ha="center", va="center", fontsize=8.5, color=CRIMSON)
b.plot([7.95, 9.45], [7.62, 7.62], color=tint(0.45), lw=2.4)
b.text(8.85, 7.25, "lncRNA, pre-mRNA", ha="center", va="center", fontsize=7.8, color=TEXT2)
b.text(8.85, 6.75, "everything else stays", ha="center", va="center", fontsize=7.8, color=TEXT2)

header(a, "Bulk RNA-seq · library type", "Poly(A) selection or rRNA depletion: choose from the question",
       dy=40, x1=3.0)
fig.savefig(OUT + "library_selection.png"); print("ok h1 library")
