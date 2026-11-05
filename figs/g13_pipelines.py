from gallery import *

# One shared look for the four W3 pipelines (after assets/img/w2/pipelines.png):
# panels = "Pre-processing", "① QC", "② Statistical test"; numbered steps with file badges;
# a crimson box at the end says what comes out. Units are inches (equal aspect).
MONO = "DejaVu Sans Mono"
BG = "#F8F9FA"
ARROW = "#9AA1AB"

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
    """items: (name, detail, badge). Returns the y below the last row."""
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
    """QC items: (name, detail) with a small crimson-tint square bullet."""
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

def three_band(name, eyebrow, title, pre, qc, test, res, res_sub,
               sub_pre, sub_qc, sub_test, qc_note=None, W=11.5):
    """Pre-processing | ① QC | ② Statistical test, side by side, then the result box."""
    gap = 0.34; w = (W - 2 * gap) / 3
    hc = [steps_h(pre), steps_h(qc) + (0.35 if qc_note else 0), steps_h(test)]
    ph = max(hc) + 1.02
    H = ph + 2.17
    fig, ax = new_fig(W, H)
    top = H - 1.05 - 0.05
    xs = [0, w + gap, 2 * (w + gap)]
    y = panel(ax, xs[0], top, w, ph, "Pre-processing", sub_pre); steps(ax, xs[0], y, w, pre)
    y = panel(ax, xs[1], top, w, ph, "① QC", sub_qc); yq = checks(ax, xs[1], y, w, qc)
    if qc_note:
        ax.text(xs[1] + 0.22, top - ph + 0.22, qc_note, ha="left", va="bottom", fontsize=8.6,
                color=TEXT2, style="italic")
    y = panel(ax, xs[2], top, w, ph, "② Statistical test", sub_test)
    steps(ax, xs[2], y, w, test, n0=len(pre) + 1)
    ym = top - 0.30
    for i in (0, 1):
        arrow(ax, (xs[i] + w + 0.03, ym), (xs[i + 1] - 0.03, ym))
    yb = top - ph
    ax.plot([xs[2] + w / 2, xs[2] + w / 2, W / 2], [yb, yb - 0.15, yb - 0.15], color=ARROW, lw=1.1)
    arrow(ax, (W / 2, yb - 0.15), (W / 2, yb - 0.27))
    result(ax, W / 2, yb - 0.27, 6.4, res, res_sub)
    header(ax, eyebrow, title)
    fig.savefig(OUT + name); print("ok", name)

# a) ATAC-seq -----------------------------------------------------------------------------------
three_band("atac_pipeline.png", "Pipeline · ATAC-seq", "ATAC-seq: from reads to open regions",
    pre=[("Sequencer output", "paired-end reads", "FASTQ"),
         ("Trim adapters", "Nextera adapters\ncutadapt / fastp", None),
         ("Align to the reference", "Bowtie2, then sort", "BAM"),
         ("Filter reads", "low MAPQ · PCR duplicates (Picard)\nchrM reads · ENCODE blacklist", None),
         ("Shift for the Tn5 offset", "+4 bp on the + strand, −5 bp on the −", "BAM")],
    qc=[("Fragment size distribution", "nucleosome-free peak + ~200 bp ladder"),
        ("TSS enrichment", "signal piles up at promoters"),
        ("FRiP", "fraction of reads in called peaks"),
        ("% chrM", "reads spent on mitochondrial DNA"),
        ("Library complexity", "unique fragments vs PCR duplicates")],
    test=[("Peak calling", "is this region more open than\nthe background? · MACS2 / MACS3", "narrowPeak"),
          ("Differential accessibility", "consensus peaks × samples counts\n→ DESeq2 / edgeR (DiffBind)", "table")],
    res="Open regions per sample — and where they differ between groups",
    res_sub="→ footprints · motif enrichment · chromatin states · GWAS overlap",
    sub_pre="reads → clean, aligned fragments", sub_qc="did the transposition work?",
    sub_test="where is the DNA open, and does it change?")

# b) ChIP-seq -----------------------------------------------------------------------------------
three_band("chip_pipeline.png", "Pipeline · ChIP-seq (also CUT&RUN)", "ChIP-seq: from reads to enriched regions",
    pre=[("Sequencer output", "ChIP and input reads", "FASTQ"),
         ("Trim and check reads", "adapter trimming · FastQC", None),
         ("Align to the reference", "Bowtie2 / BWA, then sort", "BAM"),
         ("Filter reads", "low MAPQ · PCR duplicates\nENCODE blacklist", "BAM")],
    qc=[("FRiP", "fraction of reads in called peaks"),
        ("Cross-correlation (NSC / RSC)", "strand-shifted reads line up at\nthe fragment length"),
        ("Library complexity (NRF, PBC)", "unique fragments vs PCR duplicates"),
        ("Fingerprint plot", "ChIP reads pile into few bins;\ninput reads spread evenly")],
    test=[("Peak calling", "ChIP vs input · MACS2\nnarrow mode: H3K4me3, TFs\nbroad mode: H3K27me3, H3K36me3", "narrow/broadPeak"),
          ("Differential binding", "read counts in peaks between groups\nDiffBind / csaw", "table")],
    res="Where the mark is — per sample, and where it differs between groups",
    res_sub="→ histone code · chromatin states · GWAS overlap",
    sub_pre="ChIP and input are processed identically", sub_qc="did the antibody enrich anything?",
    sub_test="enriched over input, and between groups?")

# c) Hi-C ---------------------------------------------------------------------------------------
three_band("hic_pipeline.png", "Pipeline · Hi-C", "Hi-C: from read pairs to a contact matrix",
    pre=[("Sequencer output", "paired-end reads across a ligation junction", "FASTQ"),
         ("Align each end separately", "bwa mem -SP · chimeric reads split\nat the junction", None),
         ("Parse and filter pairs", "duplicates · self-ligation · dangling ends\nsame-fragment pairs\npairtools / HiC-Pro / Juicer", "pairs"),
         ("Bin into a matrix", "contacts per pair of bins\nbin size from 5 kb to 1 Mb", ".hic/.cool"),
         ("Balance the matrix", "ICE / KR normalisation", None)],
    qc=[("Valid-pair fraction", "read pairs that are real ligation products"),
        ("Cis / trans ratio", "contacts within vs between chromosomes"),
        ("Long-range cis (> 20 kb)", "share of informative distant contacts"),
        ("Distance decay P(s)", "contact frequency falls with distance")],
    test=[("Compartments", "eigenvector of the correlation matrix → A / B", None),
          ("TADs", "insulation score → domain boundaries", None),
          ("Loops", "HiCCUPS: pixel vs its local background", None),
          ("Differential contacts", "diffHic / multiHiCcompare", "table")],
    res="A contact matrix — with compartments, TADs and loops",
    res_sub="→ which enhancer touches which promoter",
    sub_pre="read pairs → a balanced matrix", sub_qc="did proximity ligation work?",
    sub_test="structure at three scales, then differences")

# d) DNA methylation: two inputs → one β / M matrix → QC | test -------------------------------------
def methylation(name="methylation_pipeline.png", W=11.5):
    gap = 0.34; w2 = (W - gap) / 2
    bs = [("Sequencer output", "bisulfite-converted reads", "FASTQ"),
          ("Trim", "Trim Galore: adapters, low-quality ends", None),
          ("Bisulfite-aware alignment", "Bismark / bwa-meth: C→T-aware", "BAM"),
          ("Deduplicate", "WGBS only — not RRBS (reads share\nthe same cut sites)", None),
          ("Extract methylation", "methylated / unmethylated reads per CpG", "CpG report")]
    ar = [("Scanner output", "two-colour intensity per probe", "IDAT"),
          ("Detection p-value", "drop failed probes and samples", None),
          ("Background + dye-bias correction", "Noob", None),
          ("Normalise", "BMIQ / functional normalisation", None),
          ("Filter probes", "SNP-affected · cross-reactive · sex chromosomes\nminfi / SeSAMe", None)]
    qc = [("Bisulfite conversion", "non-CpG methylation, λ spike-in\narray: bisulfite control probes"),
          ("Coverage / detection p", "reads per CpG (sequencing)\ndetection p-value (array)"),
          ("Sex check", "predicted vs recorded sex"),
          ("β density and PCA", "bimodal β? do samples cluster by batch?")]
    te = [("Estimate cell-type composition", "reference-based deconvolution\n(e.g. Houseman)", None),
          ("DMP", "limma on M values · covariates: cell\ncomposition, batch, age, sex", "table"),
          ("DMR", "DMRcate / bumphunter / dmrseq", "BED")]
    p1 = max(steps_h(bs), steps_h(ar)) + 1.02
    p2 = max(steps_h(qc), steps_h(te)) + 1.02
    mh = 0.62                                   # merge box height
    H = 1.05 + 0.05 + p1 + 0.25 + mh + 0.3 + p2 + 0.25 + 0.78 + 0.08
    fig, ax = new_fig(W, H)
    top = H - 1.1
    xs = [0, w2 + gap]
    y = panel(ax, xs[0], top, w2, p1, "Pre-processing · bisulfite sequencing", "WGBS / RRBS — starts from reads")
    steps(ax, xs[0], y, w2, bs)
    y = panel(ax, xs[1], top, w2, p1, "Pre-processing · array (EPIC)", "no reads — starts from fluorescence")
    steps(ax, xs[1], y, w2, ar)
    # merge into the shared matrix
    yb = top - p1; ym = yb - 0.25
    for x in xs:
        ax.plot([x + w2 / 2, x + w2 / 2, W / 2], [yb, yb - 0.12, yb - 0.12], color=ARROW, lw=1.1)
    arrow(ax, (W / 2, yb - 0.12), (W / 2, ym))
    ax.add_patch(mp.FancyBboxPatch((W / 2 - 3.0, ym - mh), 6.0, mh, boxstyle="round,pad=0,rounding_size=0.08",
                                   fc=tint(0.12), ec=tint(0.45), lw=1, zorder=2))
    ax.text(W / 2, ym - mh / 2, "β / M value matrix  (CpG × sample)", ha="center", va="center",
            fontsize=12, fontweight="bold", color=CRIMSON, zorder=3)
    # QC | test
    t2 = ym - mh - 0.3
    arrow(ax, (W / 2, ym - mh), (W / 2, t2 + 0.02))
    y = panel(ax, xs[0], t2, w2, p2, "① QC", "did conversion and hybridisation work?"); checks(ax, xs[0], y, w2, qc)
    y = panel(ax, xs[1], t2, w2, p2, "② Statistical test", "which CpGs differ between groups?")
    steps(ax, xs[1], y, w2, te, n0=6)
    arrow(ax, (xs[0] + w2 + 0.03, t2 - 0.30), (xs[1] - 0.03, t2 - 0.30))
    yb2 = t2 - p2
    ax.plot([xs[1] + w2 / 2, xs[1] + w2 / 2, W / 2], [yb2, yb2 - 0.12, yb2 - 0.12], color=ARROW, lw=1.1)
    arrow(ax, (W / 2, yb2 - 0.12), (W / 2, yb2 - 0.25))
    result(ax, W / 2, yb2 - 0.25, 6.4, "DMPs and DMRs — which CpGs and regions differ",
           "→ annotate with chromatin states · link to genes and GWAS loci")
    header(ax, "Pipeline · DNA methylation", "DNA methylation: two inputs, one β / M matrix")
    fig.savefig(OUT + name); print("ok", name)

methylation()
