from gallery import *
import os
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "assets", "img", "w7") + "/"

# W7 pipeline, in the shared look of the W3 pipelines (g13_pipelines.py): two inputs (16S, shotgun)
# meet in one feature table → shared clean-up → ① QC | ② Statistical test → result. Units are inches.
MONO = "DejaVu Sans Mono"; BG = "#F8F9FA"; ARROW = "#9AA1AB"

def row_h(det):
    nl = det.count("\n") + 1 if det else 0
    return 0.44 + 0.21 * nl

def steps_h(items):
    return sum(row_h(d) for _, d, *_ in items)

def new_fig(W, H, top=1.05):
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0, 0, 1, (H - top) / H])
    ax.set_xlim(0, W); ax.set_ylim(0, H - top); ax.axis("off")
    return fig, ax

def panel(ax, x, ytop, w, h, title, sub=None):
    ax.add_patch(mp.FancyBboxPatch((x, ytop - h), w, h, boxstyle="round,pad=0,rounding_size=0.12",
                                   fc=BG, ec=LINE, lw=0.9, zorder=0))
    ax.text(x + 0.22, ytop - 0.30, title, ha="left", va="center", fontsize=12.5,
            fontweight="bold", color=TEXT)
    if sub:
        ax.text(x + 0.22, ytop - 0.56, sub, ha="left", va="center", fontsize=9.3, color=TEXT2)
    return ytop - (0.98 if sub else 0.72)

def steps(ax, x, y, w, items, n0=1):
    cx = x + 0.42; ys = []
    for k, (name, det, badge) in enumerate(items):
        ys.append(y)
        ax.scatter([cx], [y - 0.06], s=430, color=tint(0.16), zorder=3, linewidths=0)
        ax.text(cx, y - 0.06, str(n0 + k), ha="center", va="center", fontsize=10,
                fontweight="bold", color=CRIMSON, zorder=4)
        ax.text(cx + 0.32, y + 0.04, name, ha="left", va="center", fontsize=11,
                fontweight="semibold", color=TEXT)
        if det:
            ax.text(cx + 0.32, y - 0.09, det, ha="left", va="top", fontsize=9.2, color=TEXT2,
                    linespacing=1.35)
        if badge:
            ax.text(x + w - 0.2, y + 0.04, badge, ha="right", va="center", fontsize=8.8, color=TEXT,
                    family=MONO, zorder=4,
                    bbox=dict(boxstyle="round,pad=0.32,rounding_size=0.7", fc="white", ec=LINE, lw=0.9))
        y -= row_h(det)
    if len(ys) > 1:
        ax.plot([cx, cx], [ys[0] - 0.06, ys[-1] - 0.06], color=LINE, lw=1.1, ls=(0, (2, 2.5)), zorder=1)
    return y

def checks(ax, x, y, w, items):
    for name, det in items:
        ax.scatter([x + 0.42], [y + 0.04], marker="s", s=44, color=tint(0.55), zorder=3, linewidths=0)
        ax.text(x + 0.62, y + 0.04, name, ha="left", va="center", fontsize=10.8,
                fontweight="semibold", color=TEXT)
        if det:
            ax.text(x + 0.62, y - 0.09, det, ha="left", va="top", fontsize=9.2, color=TEXT2,
                    linespacing=1.35)
        y -= row_h(det)
    return y

def arrow(ax, p, q):
    ax.annotate("", q, p, arrowprops=dict(arrowstyle="-|>", color=ARROW, lw=1.3, mutation_scale=13,
                                          shrinkA=0, shrinkB=0), zorder=2)

def result(ax, xc, y, w, text, sub):
    ax.add_patch(mp.FancyBboxPatch((xc - w / 2, y - 0.78), w, 0.78,
                                   boxstyle="round,pad=0,rounding_size=0.06", fc=CRIMSON, ec="none", zorder=2))
    ax.text(xc, y - 0.28, text, ha="center", va="center", color="white", fontsize=12,
            fontweight="bold", zorder=3)
    ax.text(xc, y - 0.57, sub, ha="center", va="center", color=tint(0.18), fontsize=9.5, zorder=3)

W = 11.5; gap = 0.34; w2 = (W - gap) / 2
s16 = [("Sequencer output", "paired-end amplicon reads (e.g. V4)", "FASTQ"),
       ("Remove primers", "cutadapt", None),
       ("Filter and denoise", "drop low-quality reads, learn the error\nmodel → exact sequences · DADA2", None),
       ("Merge pairs, remove chimeras", "forward + reverse must overlap;\nPCR chimeras are discarded", "ASV table"),
       ("Assign taxonomy, build a tree", "classifier against SILVA / GTDB\nphylogeny for UniFrac / Faith's PD", None)]
ssg = [("Sequencer output", "reads from all DNA in the tube", "FASTQ"),
       ("Trim", "adapters, low-quality ends · fastp", None),
       ("Remove host reads", "align to the human genome, keep the rest\nKneadData / Bowtie2", None),
       ("Taxonomic profile", "marker genes (MetaPhlAn) or k-mers\n(Kraken2 + Bracken)", "taxa table"),
       ("Functional profile", "genes → pathways · HUMAnN", "pathway table")]
shared_l = [("Remove contaminants", "taxa found in blanks, or rising as\nDNA falls · decontam", None)]
shared_r = [("Handle sequencing depth", "rarefy for diversity · CLR or a model\nfor differential abundance", None)]
qc = [("Reads kept at each step", "filtered · denoised · merged · non-chimeric"),
      ("% host reads (shotgun)", "stool ≈ few %; biopsy can be > 99 %"),
      ("% classified", "reads the database could name"),
      ("Blanks and mock community", "reads in negative controls; known\ncommunity recovered as expected?"),
      ("Depth per sample", "rarefaction curves level off?")]
te = [("Alpha diversity", "richness, Shannon — within a sample", None),
      ("Beta diversity", "Bray–Curtis / UniFrac → PCoA\nPERMANOVA: report R² and p", None),
      ("Differential abundance", "which taxa differ? run ≥ 2 methods\n(e.g. ALDEx2, ANCOM-BC, MaAsLin2)", "table")]

p1 = max(steps_h(s16), steps_h(ssg)) + 1.02
mh = 0.62
sh = steps_h(shared_l) + 0.95
p2 = max(steps_h(qc), steps_h(te)) + 1.02
H = 1.05 + 0.05 + p1 + 0.25 + mh + 0.3 + sh + 0.3 + p2 + 0.25 + 0.78 + 0.08
fig, ax = new_fig(W, H)
top = H - 1.1
xs = [0, w2 + gap]
y = panel(ax, xs[0], top, w2, p1, "Pre-processing · 16S amplicon", "one marker gene — reads are near-identical copies")
steps(ax, xs[0], y, w2, s16)
y = panel(ax, xs[1], top, w2, p1, "Pre-processing · shotgun", "all DNA — microbes and host mixed")
steps(ax, xs[1], y, w2, ssg)
# merge into one feature table
yb = top - p1; ym = yb - 0.25
for x in xs:
    ax.plot([x + w2 / 2, x + w2 / 2, W / 2], [yb, yb - 0.12, yb - 0.12], color=ARROW, lw=1.1)
arrow(ax, (W / 2, yb - 0.12), (W / 2, ym))
ax.add_patch(mp.FancyBboxPatch((W / 2 - 3.3, ym - mh), 6.6, mh, boxstyle="round,pad=0,rounding_size=0.08",
                               fc=tint(0.12), ec=tint(0.45), lw=1, zorder=2))
ax.text(W / 2, ym - mh / 2, "Feature table  (taxa or pathways × samples)", ha="center", va="center",
        fontsize=12, fontweight="bold", color=CRIMSON, zorder=3)
# shared clean-up on the table
t1 = ym - mh - 0.3
arrow(ax, (W / 2, ym - mh), (W / 2, t1 + 0.02))
y = panel(ax, 0, t1, W, sh, "Pre-processing · the table", "both routes, before any statistics")
steps(ax, 0, y, w2, shared_l, n0=6)
steps(ax, xs[1], y, w2, shared_r, n0=7)
# QC | test
t2 = t1 - sh - 0.3
arrow(ax, (W / 2, t1 - sh), (W / 2, t2 + 0.02))
y = panel(ax, xs[0], t2, w2, p2, "① QC", "did each step keep the reads it should?"); checks(ax, xs[0], y, w2, qc)
y = panel(ax, xs[1], t2, w2, p2, "② Statistical test", "how diverse, how different, which taxa?")
steps(ax, xs[1], y, w2, te, n0=8)
arrow(ax, (xs[0] + w2 + 0.03, t2 - 0.30), (xs[1] - 0.03, t2 - 0.30))
yb2 = t2 - p2
ax.plot([xs[1] + w2 / 2, xs[1] + w2 / 2, W / 2], [yb2, yb2 - 0.12, yb2 - 0.12], color=ARROW, lw=1.1)
arrow(ax, (W / 2, yb2 - 0.12), (W / 2, yb2 - 0.25))
result(ax, W / 2, yb2 - 0.25, 7.0, "Which taxa and functions differ — and how much the grouping explains",
       "→ functional interpretation · artifact, confounding and reverse causation checks")
header(ax, "Pipeline · metagenomics", "Metagenomics: two routes, one feature table")
fig.savefig(OUT + "metagenomics_pipeline.png"); print("ok q1 pipeline")
