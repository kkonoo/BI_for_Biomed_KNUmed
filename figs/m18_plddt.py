from gallery import *
rng = np.random.default_rng(18)

# Reading an AlphaFold prediction: per-residue confidence (pLDDT) along an illustrative 420-residue protein with
# two folded domains joined by a disordered linker, plus a disordered tail
L = 420
x = np.arange(1, L + 1)
prof = np.full(L, 35.0)
for a_, b_ in ((20, 150), (215, 360)):                    # folded domains
    prof[a_:b_] = 93
prof[150:215] = 42; prof[360:] = 30                         # linker and tail
prof[:20] = 55
k = np.ones(9) / 9
pl = np.clip(np.convolve(prof, k, mode="same") + rng.normal(0, 2.5, L), 20, 98)
fig, ax = plt.subplots(figsize=(9.8, 3.6))
fig.subplots_adjust(left=0.07, right=0.8, top=0.72, bottom=0.17)
bands = [(90, 100, "#2B2F36", "very high (> 90)"), (70, 90, "#5E6670", "confident (70–90)"),
         (50, 70, tint(0.45), "low (50–70)"), (0, 50, CRIMSON, "very low (< 50)")]
for lo, hi, c, lab in bands:
    ax.axhspan(lo, hi, color=c, alpha=0.07, lw=0)
    m = (pl >= lo) & (pl < hi if hi < 100 else pl <= hi)
    ax.scatter(x[m], pl[m], s=5, color=c, lw=0)
    ax.text(L + 12, (lo + hi) / 2, lab, ha="left", va="center", fontsize=8.6, color=c if c != "#2B2F36" else TEXT)
ax.set_xlim(0, L); ax.set_ylim(0, 100)
ax.set_xlabel("residue", color=TEXT2, fontsize=9); ax.set_ylabel("pLDDT", color=TEXT2, fontsize=9)
for sp in ("left", "bottom"): ax.spines[sp].set_color(LINE)
ax.tick_params(colors=TEXT2, labelsize=8.3)
for a_, b_, lab in ((20, 150, "domain 1"), (215, 360, "domain 2")):
    ax.text((a_ + b_) / 2, 98.5, lab, ha="center", va="bottom", fontsize=9, color=TEXT, fontweight="semibold")
ax.text(182, 60, "linker", ha="center", va="center", fontsize=9, color=CRIMSON, fontweight="semibold")
ax.text(390, 52, "tail", ha="center", va="center", fontsize=9, color=CRIMSON, fontweight="semibold")
ax.text(0, -0.3, "High pLDDT: the fold is usually right.  Below 50: usually disordered — the drawn shape is not a "
        "structure.  How the two domains sit relative to each other is a separate score (PAE).  Illustrative profile.",
        transform=ax.transAxes, fontsize=8.4, color=TEXT2)
header(ax, "Proteomics · predicted structure", "A predicted structure is a model — read its confidence first",
       dy=14, x1=1.25)
fig.savefig(OUT + "plddt.png"); print("ok m18 plddt")
