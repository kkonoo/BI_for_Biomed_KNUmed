from gallery import *
rng = np.random.default_rng(28)

# DMP vs DMR: one CpG that differs vs a run of neighbouring CpGs that shift together
fig, axes = plt.subplots(1, 2, figsize=(9.8, 5.0), gridspec_kw={"wspace": 0.16})
DARK = "#5E6670"; GREYF = "#9AA1AB"
pos = np.array([0.0, 1.1, 1.6, 2.9, 3.5, 4.0, 4.9, 5.6, 6.1, 7.0, 8.2, 9.1])  # CpG positions (×100 bp)
base = np.array([0.82, 0.78, 0.7, 0.55, 0.42, 0.35, 0.33, 0.38, 0.45, 0.6, 0.74, 0.8])
N = 6

def draw(ax, case_mean, hl, title, sub, lab):
    lo, hi = hl
    ax.axvspan(pos[lo] - 0.3, pos[hi] + 0.3, color=FAINT, zorder=0, lw=0)
    for g, mean, col, dx in (("control", base, GREYF, -0.09), ("case", case_mean, CRIMSON, 0.09)):
        for p, m in zip(pos, mean):
            ys = np.clip(m + rng.normal(0, 0.035, N), 0, 1)
            ax.scatter(np.full(N, p + dx) + rng.uniform(-0.04, 0.04, N), ys, s=11, color=col, lw=0,
                       alpha=0.75, zorder=2)
        ax.plot(pos + dx, mean, color=col, lw=1.6, zorder=3)
    # CpG lollipops on the DNA line
    ax.plot([-0.4, 9.5], [-0.12, -0.12], color=TEXT2, lw=1.1, clip_on=False)
    for p in pos:
        ax.plot([p, p], [-0.12, -0.075], color=TEXT2, lw=0.9, clip_on=False)
        ax.plot([p], [-0.065], "o", ms=4.5, color=TEXT2, clip_on=False)
    ax.text(9.5, -0.2, "CpGs along ~1 kb", ha="right", va="top", fontsize=8.8, color=TEXT2)
    ax.text((pos[lo] + pos[hi]) / 2, 1.06, lab, ha="center", va="bottom", fontsize=10.5, fontweight="bold",
            color=CRIMSON)
    ax.set_xlim(-0.4, 9.5); ax.set_ylim(0, 1.12); ax.set_xticks([])
    ax.set_yticks([0, 0.5, 1]); ax.set_ylabel("β", color=TEXT2)
    for sp in ("left",): ax.spines[sp].set_color(LINE)
    ax.spines["bottom"].set_visible(False)
    ax.tick_params(colors=TEXT2, labelsize=8.5)
    ax.set_title(title, loc="left", fontsize=12.5, fontweight="semibold", color=TEXT, pad=26)
    ax.text(0, 1.03, sub, transform=ax.transAxes, ha="left", va="bottom", fontsize=9.3, color=TEXT2)

dmp = base.copy(); dmp[6] = 0.72
draw(axes[0], dmp, (6, 6), "DMP · one position", "a single CpG differs between groups", "DMP")
dmr = base.copy(); dmr[3:10] = np.clip(dmr[3:10] + np.array([0.12, 0.22, 0.27, 0.28, 0.26, 0.2, 0.1]), 0, 1)
draw(axes[1], dmr, (3, 9), "DMR · a region", "neighbouring CpGs shift together", "DMR")

axes[0].text(-0.4, -0.3, "tested one CpG at a time (limma on M values)\n"
             "⚠ one hit while its neighbours do not move\n   → check the probe (SNP, cross-reactive)",
             ha="left", va="top", fontsize=9, color=TEXT2, linespacing=1.45, transform=axes[0].transData)
axes[1].text(-0.4, -0.3, "DMRcate, bumphunter, dmrseq join neighbours\n"
             "✓ a run of consistent CpGs is hard to get by\n   chance or from one bad probe → more reliable",
             ha="left", va="top", fontsize=9, color=TEXT2, linespacing=1.45, transform=axes[1].transData)
# legend
axes[1].plot([6.3], [0.1], "o", ms=5, color=GREYF)
axes[1].text(6.5, 0.1, "control", ha="left", va="center", fontsize=9, color=TEXT2)
axes[1].plot([7.9], [0.1], "o", ms=5, color=CRIMSON)
axes[1].text(8.1, 0.1, "case", ha="left", va="center", fontsize=9, color=CRIMSON)

header(axes[0], "DNA methylation · differential methylation",
       "DMP = one CpG; DMR = a stretch of CpGs that moves together", dy=44, x1=2.18)
fig.savefig(OUT + "dmp_dmr.png"); print("ok g28 dmp dmr")
