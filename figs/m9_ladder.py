from gallery import *
rng = np.random.default_rng(9)

# Reading an MS2 spectrum: b and y ion ladders of the tryptic peptide LVNEVTEFAK (a serum albumin peptide).
# m/z values are computed from monoisotopic residue masses (singly charged fragments); intensities are illustrative.
RES = {"G": 57.02146, "A": 71.03711, "S": 87.03203, "P": 97.05276, "V": 99.06841, "T": 101.04768, "L": 113.08406,
       "N": 114.04293, "D": 115.02694, "K": 128.09496, "E": 129.04259, "F": 147.06841, "R": 156.10111}
PROTON, WATER = 1.00728, 18.01056
pep = "LVNEVTEFAK"
L = len(pep)
b_ion = {i: sum(RES[a] for a in pep[:i]) + PROTON for i in range(1, L)}
y_ion = {i: sum(RES[a] for a in pep[L - i:]) + WATER + PROTON for i in range(1, L)}
prec2 = (sum(RES[a] for a in pep) + WATER + 2 * PROTON) / 2

fig = plt.figure(figsize=(9.8, 5.6))
gs = fig.add_gridspec(2, 1, height_ratios=[0.42, 1], hspace=0.12, left=0.06, right=0.97, top=0.84, bottom=0.1)
DARK = "#5E6670"; GREYF = "#B8BEC7"

# --- the peptide and one break ---
t = fig.add_subplot(gs[0]); t.set_xlim(0, 10); t.set_ylim(0, 2.2); t.axis("off")
x0, dx = 1.7, 0.62
for i, a in enumerate(pep):
    x = x0 + i * dx
    t.add_patch(mp.FancyBboxPatch((x - 0.24, 0.95), 0.48, 0.5, boxstyle="round,pad=0,rounding_size=0.06",
                                  fc=DARK if i < 4 else tint(0.75), ec="none"))
    t.text(x, 1.2, a, ha="center", va="center", fontsize=11, color="white", fontweight="bold")
xb = x0 + 3.5 * dx
t.plot([xb, xb], [0.8, 1.6], color=TEXT, lw=1.2, ls=(0, (2, 1.5)))
t.text(xb, 1.75, "backbone breaks here", ha="center", va="bottom", fontsize=8.6, color=TEXT2)
t.annotate("", (x0 - 0.3, 0.68), (xb - 0.05, 0.68), arrowprops=dict(arrowstyle="-", color=DARK, lw=1.4))
t.annotate("", (xb + 0.05, 0.68), (x0 + (L - 1) * dx + 0.3, 0.68), arrowprops=dict(arrowstyle="-", color=CRIMSON, lw=1.4))
t.text((x0 - 0.3 + xb) / 2, 0.42, f"b4 = LVNE   (N-terminal piece, m/z {b_ion[4]:.1f})", ha="center", va="center",
       fontsize=8.8, color=DARK, fontweight="semibold")
t.text((xb + x0 + (L - 1) * dx + 0.3) / 2, 0.42, f"y6 = VTEFAK   (C-terminal piece, m/z {y_ion[6]:.1f})", ha="center",
       va="center", fontsize=8.8, color=CRIMSON, fontweight="semibold")
t.text(x0 + (L - 1) * dx + 0.75, 1.2, f"precursor\nm/z {prec2:.2f} (2+)", ha="left", va="center", fontsize=8.6,
       color=TEXT2, linespacing=1.2)
t.text(0, 1.2, "peptide", ha="left", va="center", fontsize=9.5, color=TEXT, fontweight="semibold")

# --- the spectrum ---
s = fig.add_subplot(gs[1])
noise = rng.uniform(120, 1080, 160)
s.vlines(noise, 0, rng.gamma(1.2, 2.5, 160), color=FAINT, lw=1, zorder=0)
yh = {i: h for i, h in zip(range(1, L), (18, 30, 62, 55, 70, 100, 64, 52, 35))}
bh_ = {i: h for i, h in zip(range(1, L), (0, 34, 22, 26, 14, 12, 9, 6, 4))}
for i in range(2, L):
    s.vlines(b_ion[i], 0, bh_[i], color=DARK, lw=1.8)
    off = -6 if i in (2, 8) else 0
    s.text(b_ion[i] + off, bh_[i] + 2, f"b{i}", ha="right" if off else "center", va="bottom", fontsize=8, color=DARK)
for i in range(1, L):
    s.vlines(y_ion[i], 0, yh[i], color=CRIMSON, lw=2)
    s.text(y_ion[i], yh[i] + 2, f"y{i}", ha="center", va="bottom", fontsize=8, color=CRIMSON, fontweight="semibold")
# read the y ladder: each gap is one residue
for i in range(3, L - 1):
    a_, b_ = y_ion[i], y_ion[i + 1]
    yy = 122
    s.annotate("", (b_, yy), (a_, yy), arrowprops=dict(arrowstyle="<->", color=CRIMSON, lw=0.8, shrinkA=0, shrinkB=0))
    s.text((a_ + b_) / 2, yy + 3, f"{pep[L - i - 1]}", ha="center", va="bottom", fontsize=10.5, color=CRIMSON,
           fontweight="bold")
    s.text((a_ + b_) / 2, yy - 4, f"{b_ - a_:.2f}", ha="center", va="top", fontsize=7, color=TEXT2)
s.set_xlim(100, 1100); s.set_ylim(0, 140)
s.set_yticks([]); s.set_xlabel("m/z", color=TEXT2, fontsize=9)
for sp in ("left", "top", "right"): s.spines[sp].set_visible(False)
s.spines["bottom"].set_color(LINE); s.tick_params(colors=TEXT2, labelsize=8.3)
s.text(100, 134, "gaps between y ions = one amino acid each →\nread right to left: …V N E V T E…", ha="left",
       va="top", fontsize=8.8, color=TEXT, linespacing=1.25)
s.text(1100, 80, "red = y ions (C-terminal pieces)\ngrey = b ions (N-terminal pieces)\nintensities illustrative",
       ha="right", va="top", fontsize=8.4, color=TEXT2, linespacing=1.3)
header(t, "Proteomics · MS2 spectrum", "The sequence is in the gaps between fragment peaks", dy=8, x1=1)
fig.savefig(OUT + "ms2_ladder.png")
print("ok m9 ladder", round(prec2, 3), {k: round(v, 2) for k, v in y_ion.items()})
