from gallery import *
rng = np.random.default_rng(27)

# DNA methylation QC on one page: conversion, signal per CpG, probe filters, sample checks
fig = plt.figure(figsize=(9.8, 9.4))
gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.14, left=0.01, right=0.99, top=0.86, bottom=0.04)
MONO = "DejaVu Sans Mono"; DARK = "#5E6670"; GREYF = "#B8BEC7"

def panel(ax, t, sub):
    ax.set_title(t, loc="left", fontsize=12.5, fontweight="semibold", color=TEXT, pad=22)
    ax.text(0, 1.02, sub, transform=ax.transAxes, ha="left", va="bottom", fontsize=9.3, color=TEXT2)
def style(ax):
    for sp in ("left", "bottom"): ax.spines[sp].set_color(LINE)
    ax.tick_params(colors=TEXT2, labelsize=8)

# (a) bisulfite conversion efficiency ----------------------------------------------------------------------
a = fig.add_subplot(gs[0, 0]); a.set_xlim(0, 10); a.set_ylim(-0.4, 4.4); a.axis("off")
panel(a, "(a) Bisulfite conversion", "C outside CpG is almost never methylated → it must read T")
ref = "TCAGCGTCCACGTC"
cpg = {4, 10}                                   # C of each CpG (0-based)
def row(y, lab, seq, cols, lc=TEXT):
    a.text(0.0, y, lab, ha="left", va="center", fontsize=9.3, fontweight="semibold", color=lc)
    for i, (ch, cc) in enumerate(zip(seq, cols)):
        a.text(2.25 + i * 0.5, y, ch, ha="center", va="center", fontsize=11.5, family=MONO, color=cc,
               fontweight="bold" if cc != TEXT2 else "normal")
for i in cpg:
    a.add_patch(mp.FancyBboxPatch((2.25 + i * 0.5 - 0.22, 0.55), 0.94, 3.25, boxstyle="round,pad=0,rounding_size=0.08",
                                  fc=FAINT, ec="none", zorder=0))
row(3.45, "reference", ref, [TEXT2] * len(ref))
good = "".join("C" if (ch == "C" and i in cpg) else ("T" if ch == "C" else ch) for i, ch in enumerate(ref))
row(2.45, "good read", good, [CRIMSON if i in cpg else (TEXT if ref[i] == "C" else TEXT2) for i in range(len(ref))])
bad = list(good); bad[1] = "C"; bad[8] = "C"; bad = "".join(bad)
row(1.45, "failed read", bad, [CRIMSON if i in cpg else (DARK if i in (1, 8) else (TEXT if ref[i] == "C" else TEXT2))
                               for i in range(len(ref))], lc=TEXT2)
for i in (1, 8):
    a.add_patch(mp.Circle((2.25 + i * 0.5, 1.45), 0.2, fc="none", ec=DARK, lw=1.2))
a.text(2.25 + 4.5 * 0.5, 3.95, "CpG", ha="center", va="center", fontsize=8.8, color=TEXT2)
a.text(2.25 + 10.5 * 0.5, 3.95, "CpG", ha="center", va="center", fontsize=8.8, color=TEXT2)
a.text(0.0, 0.55, "unconverted C outside CpG → looks like methylation everywhere", ha="left", va="center",
       fontsize=8.8, color=TEXT2)
a.text(0.0, 0.0, "conversion rate = non-CpG C read as T / all non-CpG C  →  should be ≥ 99%",
       ha="left", va="center", fontsize=9.3, color=CRIMSON, fontweight="semibold")

# (b) enough signal per CpG -------------------------------------------------------------------------------
b = fig.add_subplot(gs[0, 1]); b.axis("off")
panel(b, "(b) Enough signal per CpG", "sequencing: coverage · array: detection p-value")
b1 = b.inset_axes([0.0, 0.56, 1.0, 0.36])
b1.set_xlim(-0.02, 1.25); b1.set_ylim(-0.6, 2.6); b1.axis("off")
b1.plot([0, 1], [-0.25, -0.25], color=LINE, lw=1)
for v in (0, 0.5, 1): b1.text(v, -0.45, f"{v:g}", ha="center", va="top", fontsize=8, color=TEXT2)
b1.text(0.5, 2.55, "β of one CpG (true β = 0.6)", ha="center", va="top", fontsize=8.8, color=TEXT2)
for y, n, k, lab in ((1.5, 3, 2, "3 reads"), (0.45, 30, 18, "30 reads")):
    p = k / n; se = 1.96 * np.sqrt(p * (1 - p) / n)
    lo, hi = max(0, p - se), min(1, p + se)
    if n == 3: lo, hi = 0.09, 0.99              # exact (Clopper–Pearson) interval for 2/3
    good = n == 30
    b1.plot([lo, hi], [y, y], color=CRIMSON if good else GREYF, lw=3, solid_capstyle="round")
    b1.plot([p], [y], "o", ms=7, color=CRIMSON if good else DARK, mec="white", mew=1.2, zorder=3)
    b1.text(1.05, y, lab + ("  ✓" if good else "  ✕"), ha="left", va="center", fontsize=9, fontweight="semibold",
            color=CRIMSON if good else DARK)
b2 = b.inset_axes([0.08, 0.0, 0.88, 0.38]); style(b2)
xs = np.linspace(0, 10, 400)
neg = np.exp(-0.5 * ((xs - 2.0) / 0.7) ** 2)
b2.fill_between(xs, neg, color=GREYF, lw=0)
b2.text(0.15, 1.45, "negative-control\nprobes (background)", ha="left", va="top", fontsize=8.5, color=TEXT2,
        linespacing=1.2)
for xv, ok in ((7.0, True), (3.0, False)):
    b2.axvline(xv, color=CRIMSON if ok else DARK, lw=1.6, ls="-" if ok else (0, (3, 2)))
b2.text(7.15, 0.75, "probe A: far above\n→ p < 0.01 ✓", ha="left", va="center", fontsize=8.6, color=CRIMSON,
        fontweight="semibold", linespacing=1.25)
b2.text(3.2, 0.45, "probe B: inside\nbackground → drop", ha="left", va="center", fontsize=8.6, color=DARK,
        linespacing=1.25)
b2.set_xlim(0, 10); b2.set_ylim(0, 1.5); b2.set_yticks([]); b2.set_xticks([])
b2.set_xlabel("total signal (array)", color=TEXT2, fontsize=8.8)

# (c) probe filters (array) -------------------------------------------------------------------------------
c = fig.add_subplot(gs[1, 0]); c.set_xlim(0, 10); c.set_ylim(-0.5, 4.4); c.axis("off")
panel(c, "(c) Probe filters (array)", "drop probes that measure something other than methylation")
# SNP probe: beta splits by genotype
c1 = c.inset_axes([0.0, 0.66, 0.3, 0.3]); style(c1)
for k, (mu, n) in enumerate(((0.03, 14), (0.5, 12), (0.96, 6))):
    c1.scatter(rng.uniform(0.1, 0.9, n), np.clip(rng.normal(mu, 0.025, n), 0, 1), s=9, color=DARK, lw=0)
c1.set_ylim(-0.08, 1.08); c1.set_xlim(0, 1); c1.set_xticks([]); c1.set_yticks([0, 0.5, 1])
c1.set_ylabel("β", color=TEXT2, fontsize=8.5)
c.text(3.6, 3.95, "SNP under the probe", ha="left", va="center", fontsize=10, fontweight="semibold", color=TEXT)
c.text(3.6, 3.45, "β clusters at 0 / 0.5 / 1 = the genotype,\nnot methylation", ha="left", va="center",
       fontsize=8.8, color=TEXT2, linespacing=1.3)
# cross-reactive
y = 1.85
c.plot([0.3, 2.6], [y + 0.25, y + 0.25], color=TEXT2, lw=3)
c.text(1.45, y + 0.55, "one probe", ha="center", va="center", fontsize=8.5, color=TEXT2)
for k, yy in enumerate((y - 0.35, y - 0.7, y - 1.05)):
    c.plot([0.4 + k * 0.7, 1.0 + k * 0.7], [yy, yy], color=CRIMSON if k == 0 else GREYF, lw=3)
    c.annotate("", (0.7 + k * 0.7, yy + 0.08), (1.45, y + 0.15), arrowprops=dict(arrowstyle="-|>", color=LINE,
                                                                                lw=0.9, mutation_scale=7))
c.text(3.6, y + 0.15, "cross-reactive probe", ha="left", va="center", fontsize=10, fontweight="semibold",
       color=TEXT)
c.text(3.6, y - 0.35, "binds several places in the genome →\nsignal is a mixture", ha="left", va="center",
       fontsize=8.8, color=TEXT2, linespacing=1.3)
c.text(3.6, 0.2, "sex-chromosome probes", ha="left", va="center", fontsize=10, fontweight="semibold", color=TEXT)
c.text(3.6, -0.25, "X inactivation in women → differ by sex", ha="left", va="center", fontsize=8.8, color=TEXT2)
c.text(0.3, 0.2, "chrX  chrY", ha="left", va="center", fontsize=10, family=MONO, color=DARK, fontweight="bold")

# (d) sample checks ------------------------------------------------------------------------------------------
d = fig.add_subplot(gs[1, 1]); d.axis("off")
panel(d, "(d) Sample checks", "catch swapped, failed, and batch-shifted samples")
w = 0.3
d1 = d.inset_axes([0.0, 0.3, w, 0.62]); style(d1)
d1.scatter(rng.normal(0.25, 0.03, 12), rng.normal(0.9, 0.05, 12), s=12, color=DARK, lw=0)    # male
d1.scatter(rng.normal(0.45, 0.03, 12), rng.normal(0.1, 0.05, 12), s=12, color=tint(0.45), lw=0)  # female
d1.scatter([0.25], [0.88], s=40, facecolor="none", edgecolor=CRIMSON, lw=1.6)
d1.text(0.27, 0.7, "labelled F\nbut looks M", ha="left", va="top", fontsize=7.8, color=CRIMSON, linespacing=1.2)
d1.set_xlim(0.15, 0.55); d1.set_ylim(-0.1, 1.15); d1.set_xticks([]); d1.set_yticks([])
d1.set_xlabel("chrX mean β", color=TEXT2, fontsize=8.5); d1.set_ylabel("chrY signal", color=TEXT2, fontsize=8.5)
d1.set_title("sex check", fontsize=9.5, color=TEXT, fontweight="semibold", loc="left", pad=4)
d2 = d.inset_axes([0.36, 0.3, w, 0.62]); style(d2)
bx = np.linspace(0, 1, 300)
for k in range(6):
    yv = 1.6 * np.exp(-0.5 * ((bx - 0.07) / 0.05) ** 2) + 1.3 * np.exp(-0.5 * ((bx - 0.88) / 0.07) ** 2) \
         + 0.12 + 0.05 * k * 0
    d2.plot(bx, yv * (1 + 0.05 * rng.standard_normal()), color=GREYF, lw=1)
d2.plot(bx, 0.55 + 0.5 * np.exp(-0.5 * ((bx - 0.5) / 0.22) ** 2), color=CRIMSON, lw=1.6)
d2.text(0.5, 1.15, "outlier", ha="center", va="bottom", fontsize=8, color=CRIMSON)
d2.set_xticks([0, 1]); d2.set_yticks([]); d2.set_xlabel("β", color=TEXT2, fontsize=8.5)
d2.set_title("β distribution", fontsize=9.5, color=TEXT, fontweight="semibold", loc="left", pad=4)
d3 = d.inset_axes([0.72, 0.3, w - 0.02, 0.62]); style(d3)
d3.scatter(rng.normal(-1, 0.35, 10), rng.normal(0, 0.6, 10), s=12, color=DARK, lw=0)
d3.scatter(rng.normal(1.2, 0.35, 10), rng.normal(0, 0.6, 10), s=12, color=tint(0.55), lw=0)
d3.text(-1, 1.45, "batch 1", ha="center", va="bottom", fontsize=8, color=DARK)
d3.text(1.2, 1.45, "batch 2", ha="center", va="bottom", fontsize=8, color=CRIMSON)
d3.set_xlim(-2.2, 2.4); d3.set_ylim(-1.8, 2.0); d3.set_xticks([]); d3.set_yticks([])
d3.set_xlabel("PC1", color=TEXT2, fontsize=8.5); d3.set_ylabel("PC2", color=TEXT2, fontsize=8.5)
d3.set_title("PCA", fontsize=9.5, color=TEXT, fontweight="semibold", loc="left", pad=4)
d.text(0.0, 0.08, "sample swap · failed sample · batch → fix or drop before testing", transform=d.transAxes,
       ha="left", va="center", fontsize=9, color=TEXT2)

header(a, "DNA methylation QC", "Check the conversion, the signal, the probes, then the samples", dy=40, x1=2.03)
fig.savefig(OUT + "methylation_qc.png"); print("ok g27 methylation qc")
