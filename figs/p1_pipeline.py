from gallery import *
import os

# W6 proteomics pipeline in the same look as the W3 pipelines (g13_pipelines.py):
# panels "Pre-processing", "① QC", "② Statistical test"; numbered steps with file badges;
# a crimson box at the end says what comes out. Units are inches.
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "assets", "img", "w6") + "/"
MONO = "DejaVu Sans Mono"; BG = "#F8F9FA"; ARROW = "#9AA1AB"

def row_h(det):
    nl = det.count("\n") + 1 if det else 0
    return 0.44 + 0.21 * nl

def steps_h(items):
    return sum(row_h(d) for _, d, *_ in items)

def panel(ax, x, ytop, w, h, title, sub):
    ax.add_patch(mp.FancyBboxPatch((x, ytop - h), w, h, boxstyle="round,pad=0,rounding_size=0.12",
                                   fc=BG, ec=LINE, lw=0.9, zorder=0))
    ax.text(x + 0.22, ytop - 0.30, title, ha="left", va="center", fontsize=12.5, fontweight="bold", color=TEXT)
    ax.text(x + 0.22, ytop - 0.56, sub, ha="left", va="center", fontsize=9.3, color=TEXT2)
    return ytop - 0.98

def steps(ax, x, y, w, items, n0=1):
    cx = x + 0.42; ys = []
    for k, (name, det, badge) in enumerate(items):
        ys.append(y)
        ax.scatter([cx], [y - 0.06], s=430, color=tint(0.16), zorder=3, linewidths=0)
        ax.text(cx, y - 0.06, str(n0 + k), ha="center", va="center", fontsize=10, fontweight="bold",
                color=CRIMSON, zorder=4)
        ax.text(cx + 0.32, y + 0.04, name, ha="left", va="center", fontsize=11, fontweight="semibold", color=TEXT)
        if det:
            ax.text(cx + 0.32, y - 0.09, det, ha="left", va="top", fontsize=9.2, color=TEXT2, linespacing=1.35)
        if badge:
            ax.text(x + w - 0.2, y + 0.04, badge, ha="right", va="center", fontsize=8.8, color=TEXT, family=MONO,
                    zorder=4, bbox=dict(boxstyle="round,pad=0.32,rounding_size=0.7", fc="white", ec=LINE, lw=0.9))
        y -= row_h(det)
    if len(ys) > 1:
        ax.plot([cx, cx], [ys[0] - 0.06, ys[-1] - 0.06], color=LINE, lw=1.1, ls=(0, (2, 2.5)), zorder=1)

def checks(ax, x, y, items):
    for name, det in items:
        ax.scatter([x + 0.42], [y + 0.04], marker="s", s=44, color=tint(0.55), zorder=3, linewidths=0)
        ax.text(x + 0.62, y + 0.04, name, ha="left", va="center", fontsize=10.8, fontweight="semibold", color=TEXT)
        if det:
            ax.text(x + 0.62, y - 0.09, det, ha="left", va="top", fontsize=9.2, color=TEXT2, linespacing=1.35)
        y -= row_h(det)

def arrow(ax, p, q):
    ax.annotate("", q, p, arrowprops=dict(arrowstyle="-|>", color=ARROW, lw=1.3, mutation_scale=13,
                                          shrinkA=0, shrinkB=0), zorder=2)

pre = [("Instrument output", "one LC–MS/MS run per sample\n(per set of samples for TMT)", ".raw"),
       ("Database search", "match each MS2 spectrum to a peptide\nMaxQuant · FragPipe (DDA) · DIA-NN (DIA)", "PSMs"),
       ("Target–decoy FDR", "keep matches at 1% FDR — for PSMs,\npeptides and proteins", None),
       ("Protein inference", "peptides → protein groups", "groups"),
       ("Quantify proteins", "peptide intensities → one value per protein\nMaxLFQ (label-free) · reporter ions (TMT)", None),
       ("Log₂, normalise, filter", "median centring · keep proteins seen in\n≥ 70% of a group · impute only if needed", "matrix")]
qc = [("LC–MS run", "total ion current · retention-time\ndrift · mass error"),
      ("Identification", "PSMs, peptides and proteins per run ·\nmissed cleavages (= digestion)"),
      ("Contaminants", "keratins, trypsin; blood proteins\nin a tissue sample"),
      ("Quantification", "% missing per sample · CV of a pooled\nQC sample injected again and again"),
      ("Sample overview", "intensity distributions · correlation · PCA")]
test = [("Differential abundance", "log₂ intensity ~ group (+ batch)\nlimma · MSstats", "table"),
        ("Multiple testing", "BH FDR over 3,000–9,000 proteins", None)]

W = 11.5; gap = 0.34; w = (W - 2 * gap) / 3
ph = max(steps_h(pre), steps_h(qc), steps_h(test)) + 1.02
H = ph + 2.17
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, (H - 1.05) / H]); ax.set_xlim(0, W); ax.set_ylim(0, H - 1.05); ax.axis("off")
top = H - 1.05 - 0.05
xs = [0, w + gap, 2 * (w + gap)]
steps(ax, xs[0], panel(ax, xs[0], top, w, ph, "Pre-processing", "spectra → a protein × sample matrix"), w, pre)
checks(ax, xs[1], panel(ax, xs[1], top, w, ph, "① QC", "did the runs and the search work?"), qc)
steps(ax, xs[2], panel(ax, xs[2], top, w, ph, "② Statistical test", "which proteins differ between groups?"), w,
      test, n0=len(pre) + 1)
for i in (0, 1):
    arrow(ax, (xs[i] + w + 0.03, top - 0.30), (xs[i + 1] - 0.03, top - 0.30))
yb = top - ph
ax.plot([xs[2] + w / 2, xs[2] + w / 2, W / 2], [yb, yb - 0.15, yb - 0.15], color=ARROW, lw=1.1)
arrow(ax, (W / 2, yb - 0.15), (W / 2, yb - 0.27))
y = yb - 0.27; bw = 6.6
ax.add_patch(mp.FancyBboxPatch((W / 2 - bw / 2, y - 0.78), bw, 0.78, boxstyle="round,pad=0,rounding_size=0.06",
                               fc=CRIMSON, ec="none", zorder=2))
ax.text(W / 2, y - 0.28, "Proteins that differ — each with its protein group and unique peptides", ha="center",
        va="center", color="white", fontsize=12, fontweight="bold", zorder=3)
ax.text(W / 2, y - 0.57, "→ enrichment against the detected proteome · phosphosites · pQTLs", ha="center",
        va="center", color=tint(0.18), fontsize=9.5, zorder=3)
header(ax, "Pipeline · bottom-up proteomics", "Proteomics: from spectra to proteins that differ")
fig.savefig(OUT + "proteomics_pipeline.png"); print("ok p1 pipeline")
