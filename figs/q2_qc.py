from gallery import *
import os
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "assets", "img", "w7") + "/"
rng = np.random.default_rng(31)

# QC by pre-processing step (illustrative numbers):
# (a) 16S: share of input reads kept after each DADA2 step; one library loses most reads at merging
#     (its reads were trimmed too short for the two ends to overlap).
# (b) shotgun: share of host reads in stool vs tissue biopsy, and the microbial reads that are left.
# (c) negative controls: reads per library — samples vs extraction blanks; one blank sits among the samples.
DARK = "#5E6670"; GREYF = "#B8BEC7"
fig, axes = plt.subplots(1, 3, figsize=(11.6, 4.3), gridspec_kw={"wspace": 0.6, "width_ratios": [1.25, 1.05, 0.9]})

def tidy(a):
    for sp in ("top", "right"): a.spines[sp].set_visible(False)
    for sp in ("left", "bottom"): a.spines[sp].set_color(LINE)
    a.tick_params(colors=TEXT2, labelsize=8.3); a.minorticks_off()

# (a) DADA2 read tracking ------------------------------------------------------------------------
a = axes[0]
stepsx = ["input", "filtered", "denoised", "merged", "non-\nchimeric"]
x = np.arange(len(stepsx))
for k in range(7):
    f = 100 * np.cumprod([1, rng.uniform(0.85, 0.92), rng.uniform(0.96, 0.99),
                          rng.uniform(0.90, 0.95), rng.uniform(0.93, 0.97)])
    a.plot(x, f, color=GREYF, lw=1.3, marker="o", ms=3.5, zorder=2)
bad = 100 * np.cumprod([1, 0.88, 0.97, 0.22, 0.95])
a.plot(x, bad, color=CRIMSON, lw=1.8, marker="o", ms=4.5, zorder=3)
a.text(3.12, bad[3] - 2, "S5: reads too short\nto overlap → lost\nat merging", ha="left", va="top",
       fontsize=8.6, color=CRIMSON, fontweight="semibold", linespacing=1.15)
a.text(4.1, 68, "other samples:\n~70 % kept", ha="right", va="top", fontsize=8.6, color=DARK, linespacing=1.15)
a.set_xticks(x); a.set_xticklabels(stepsx, fontsize=8.6, color=TEXT)
a.set_ylim(0, 105); a.set_xlim(-0.25, 4.45)
a.set_ylabel("% of input reads kept", color=TEXT2, fontsize=9)
tidy(a)
a.set_title("(a) 16S: reads kept at each step", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
a.text(0, 1.03, "8 libraries through DADA2; grey = typical", transform=a.transAxes, fontsize=9, color=TEXT2)

# (b) host reads --------------------------------------------------------------------------------
b = axes[1]
labs = ["stool 1", "stool 2", "stool 3", "biopsy 1", "biopsy 2", "biopsy 3"]
host = np.array([0.6, 1.4, 2.1, 97.2, 99.1, 99.6])
yy = np.arange(len(labs))[::-1]
b.barh(yy, host, color=GREYF, height=0.62, label="host (human) reads")
b.barh(yy, 100 - host, left=host, color=CRIMSON, height=0.62, label="microbial reads")
for yv, h in zip(yy, host):
    left_m = 20 * (100 - h) / 100            # million microbial reads out of 20 M
    txt = f"{left_m:.1f} M" if left_m >= 1 else f"{left_m * 1000:.0f} k"
    b.text(101.5, yv, txt, ha="left", va="center", fontsize=8.6,
           color=CRIMSON if h > 50 else TEXT, fontweight="semibold" if h > 50 else "normal")
b.set_yticks(yy); b.set_yticklabels(labs, fontsize=8.8, color=TEXT)
b.set_xlim(0, 100); b.set_xticks([0, 50, 100]); b.set_xticklabels(["0%", "50%", "100%"])
b.set_xlabel("share of reads (20 M per library)", color=TEXT2, fontsize=9)
b.legend(loc="upper center", bbox_to_anchor=(0.45, -0.17), ncol=2, frameon=False, fontsize=8.4,
         handlelength=1, columnspacing=1)
tidy(b); b.tick_params(axis="y", length=0)
b.set_title("(b) Shotgun: host reads", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
b.text(0, 1.03, "number = microbial reads left", transform=b.transAxes, fontsize=9, color=TEXT2)

# (c) negative controls ---------------------------------------------------------------------------
c = axes[2]
samp = np.exp(rng.normal(np.log(45000), 0.35, 24))
blank = np.exp(rng.normal(np.log(600), 0.6, 5))
odd = 28000                                     # a blank that looks like a sample
c.scatter(rng.uniform(-0.18, 0.18, len(samp)), samp, s=18, color=DARK, lw=0, zorder=3)
c.scatter(1 + rng.uniform(-0.12, 0.12, len(blank)), blank, s=22, color=GREYF, lw=0, zorder=3)
c.scatter([1.0], [odd], s=34, color=CRIMSON, lw=0, zorder=4)
c.text(1.2, odd, "a blank with as many\nreads as a sample →\ncontamination or\nwell-to-well leakage",
       ha="left", va="center", fontsize=8.4, color=CRIMSON, fontweight="semibold", linespacing=1.15)
c.text(1.2, np.median(blank), "typical blank:\nfew reads", ha="left", va="center", fontsize=8.4, color=DARK,
       linespacing=1.1)
c.set_yscale("log"); c.set_ylim(80, 3e5)
c.set_yticks([1e2, 1e3, 1e4, 1e5]); c.set_yticklabels(["100", "1k", "10k", "100k"])
c.set_xticks([0, 1]); c.set_xticklabels(["samples", "extraction\nblanks"], fontsize=8.8, color=TEXT)
c.set_xlim(-0.5, 2.3)
c.set_ylabel("reads per library", color=TEXT2, fontsize=9)
tidy(c); c.tick_params(axis="x", length=0)
c.set_title("(c) Negative controls", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
c.text(0, 1.03, "blanks go through every step", transform=c.transAxes, fontsize=9, color=TEXT2)

header(a, "Metagenomics · QC", "Check each step: what it kept, what it removed, what the blanks contain",
       dy=40, x1=3.75)
fig.savefig(OUT + "meta_qc.png"); print("ok q2 qc")
