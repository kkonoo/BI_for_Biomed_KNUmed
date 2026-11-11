from gallery import *
rng = np.random.default_rng(11)

# Bottom-up workflow: what the instrument records (peptide and fragment masses) vs what software concludes (protein)
W, H = 9.8, 4.0
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, (H - 0.8) / H]); ax.set_xlim(0, W); ax.set_ylim(0, H - 0.8); ax.axis("off")
DARK = "#5E6670"; GREYF = "#B8BEC7"
N = 7; CW, CH, GAP = 1.2, 1.3, 0.18
X0 = (W - (N * CW + (N - 1) * GAP)) / 2
YC = 2.08
steps = [("Proteins", "a complex mixture"), ("Digest", "trypsin cuts after K/R"), ("LC", "peptides elute over time"),
         ("MS1", "peptide masses;\none is picked"), ("MS2", "masses of its\nfragments"),
         ("Search", "match to database\npeptides"), ("Infer", "which protein\nwas there?")]
xs = [X0 + i * (CW + GAP) for i in range(N)]

def frame(i, inferred=False):
    x0 = xs[i]
    ax.add_patch(mp.FancyBboxPatch((x0, YC - CH / 2), CW, CH, boxstyle="round,pad=0,rounding_size=0.08",
                                   fc="white", ec=CRIMSON if inferred else LINE, lw=1.1 if inferred else 0.9,
                                   ls=(0, (3, 2)) if inferred else "-"))
    nm, d = steps[i]
    ax.text(x0 + CW / 2, YC - CH / 2 - 0.2, f"{i + 1}  {nm}", ha="center", va="center", fontsize=10,
            fontweight="semibold", color=CRIMSON if inferred else TEXT)
    ax.text(x0 + CW / 2, YC - CH / 2 - 0.42, d, ha="center", va="top", fontsize=8, color=TEXT2, linespacing=1.2)
    if i < N - 1:
        ax.annotate("", (xs[i + 1] - 0.02, YC), (x0 + CW + 0.02, YC),
                    arrowprops=dict(arrowstyle="-|>", color=GREYF, lw=1, mutation_scale=8, shrinkA=0, shrinkB=0))
    return x0

# 1 proteins
x0 = frame(0)
for (dx, dy, r, c) in ((0.32, 0.3, 0.2, DARK), (0.78, 0.38, 0.15, GREYF), (0.55, -0.2, 0.24, GREYF),
                       (0.9, -0.3, 0.11, DARK), (0.25, -0.3, 0.1, GREYF)):
    ax.add_patch(mp.Circle((x0 + dx, YC + dy), r, fc=c, ec="none"))
# 2 digest: one protein chain cut into peptides
x0 = frame(1)
cuts = [0.12, 0.38, 0.58, 0.88, 1.08]
for k in range(len(cuts) - 1):
    ya = YC + 0.32 - 0.0
    ax.plot([x0 + cuts[k] + 0.03, x0 + cuts[k + 1] - 0.03], [YC + 0.3, YC + 0.3], color=DARK, lw=4, solid_capstyle="butt")
for k in range(4):
    a_, b_ = cuts[k], cuts[k + 1]
    yy = YC - 0.2 - 0.12 * (k % 2)
    ax.plot([x0 + a_ + 0.03, x0 + b_ - 0.03], [yy, yy], color=tint(0.55) if k == 1 else GREYF, lw=4,
            solid_capstyle="butt")
for c in cuts[1:-1]:
    ax.plot([x0 + c, x0 + c], [YC + 0.17, YC + 0.43], color=CRIMSON, lw=1)
ax.annotate("", (x0 + CW / 2, YC + 0.02), (x0 + CW / 2, YC + 0.16),
            arrowprops=dict(arrowstyle="-|>", color=GREYF, lw=0.8, mutation_scale=6))
# 3 LC chromatogram
x0 = frame(2)
t = np.linspace(0, 1, 200)
cg = sum(h * np.exp(-0.5 * ((t - c) / 0.035) ** 2) for c, h in ((0.15, 0.4), (0.3, 0.9), (0.42, 0.55), (0.55, 1),
                                                                 (0.68, 0.35), (0.8, 0.7), (0.9, 0.3)))
ax.fill_between(x0 + 0.1 + t * (CW - 0.2), YC - 0.4, YC - 0.4 + cg * 0.8, color=GREYF, lw=0)
ax.plot([x0 + 0.1, x0 + CW - 0.1], [YC - 0.4, YC - 0.4], color=LINE, lw=0.8)
ax.text(x0 + CW - 0.1, YC - 0.46, "time", ha="right", va="top", fontsize=7, color=TEXT2)
# 4 MS1: precursors, one picked
x0 = frame(3)
pk = [(0.12, 0.35), (0.25, 0.8), (0.4, 0.5), (0.52, 1.0), (0.66, 0.3), (0.78, 0.65), (0.9, 0.25)]
for c, h in pk:
    ax.plot([x0 + 0.1 + c * (CW - 0.2)] * 2, [YC - 0.4, YC - 0.4 + h * 0.8], color=CRIMSON if c == 0.52 else GREYF,
            lw=2 if c == 0.52 else 1.4, solid_capstyle="butt")
ax.plot([x0 + 0.1, x0 + CW - 0.1], [YC - 0.4, YC - 0.4], color=LINE, lw=0.8)
ax.text(x0 + CW - 0.1, YC - 0.46, "m/z", ha="right", va="top", fontsize=7, color=TEXT2)
# 5 MS2: fragment ladder
x0 = frame(4)
for k, c in enumerate(np.linspace(0.08, 0.92, 9)):
    h = rng.uniform(0.3, 1.0)
    ax.plot([x0 + 0.1 + c * (CW - 0.2)] * 2, [YC - 0.4, YC - 0.4 + h * 0.8], color=tint(0.75) if k % 2 else DARK,
            lw=1.4, solid_capstyle="butt")
ax.plot([x0 + 0.1, x0 + CW - 0.1], [YC - 0.4, YC - 0.4], color=LINE, lw=0.8)
ax.text(x0 + CW - 0.1, YC - 0.46, "m/z", ha="right", va="top", fontsize=7, color=TEXT2)
# 6 search
x0 = frame(5, True)
for k, (sq, sc) in enumerate((("LVNEVTEFAK", True), ("AEFAEVSK", False), ("YLYEIAR", False), ("…", False))):
    ax.text(x0 + 0.12, YC + 0.4 - k * 0.27, sq, ha="left", va="center", fontsize=7.8, family="DejaVu Sans Mono",
            color=CRIMSON if sc else TEXT2, fontweight="bold" if sc else "normal")
ax.text(x0 + CW - 0.1, YC + 0.4, "✓", ha="right", va="center", fontsize=9, color=CRIMSON)
# 7 infer
x0 = frame(6, True)
ax.text(x0 + CW / 2, YC + 0.18, "albumin", ha="center", va="center", fontsize=10, fontweight="bold", color=CRIMSON)
ax.text(x0 + CW / 2, YC - 0.15, "(a conclusion,\nnot an observation)", ha="center", va="center", fontsize=7.4,
        color=TEXT2, linespacing=1.15)

# bands: what is measured vs what is concluded
def band(i0, i1, y, text, col):
    a_, b_ = xs[i0], xs[i1] + CW
    ax.plot([a_, b_], [y, y], color=col, lw=1.4)
    for xx in (a_, b_): ax.plot([xx, xx], [y, y + 0.1], color=col, lw=1.4)
    ax.text((a_ + b_) / 2, y - 0.17, text, ha="center", va="center", fontsize=9, color=col, fontweight="semibold")
band(0, 1, 0.36, "sample preparation", TEXT2)
band(2, 4, 0.36, "what the instrument records: masses", DARK)
band(5, 6, 0.36, "what software concludes", CRIMSON)
header(ax, "Proteomics · bottom-up workflow", "You never observe a protein — you infer it from peptides", x=X0 / W,
       dy=0, x1=1 - X0 / W)
fig.savefig(OUT + "bottomup_msms.png"); print("ok m2 workflow")
