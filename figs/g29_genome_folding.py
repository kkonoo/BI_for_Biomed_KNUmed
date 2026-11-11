from gallery import *
from matplotlib.colors import LinearSegmentedColormap
rng = np.random.default_rng(29)

# Folding levels of the genome, from DNA to chromosome territory, and what each looks like in Hi-C
CM = LinearSegmentedColormap.from_list("hic", ["#FFFFFF", tint(0.25), tint(0.6), CRIMSON, "#5E1220"])
FW, FH = 9.8, 6.0
fig = plt.figure(figsize=(FW, FH))
DARK = "#5E6670"; PROT = "#C9CDD3"; GREYF = "#B8BEC7"
N = 6; W = 0.142; GAP = (1 - N * W) / (N - 1)
X = [k * (W + GAP) for k in range(N)]
CY, CH = 0.53, 0.25            # cartoon row
MY, MH = 0.06, W * FW / FH      # matrix row (square)

levels = [
    ("DNA", "2 nm", "double helix"),
    ("Nucleosome", "11 nm · ~200 bp", "DNA wound on histones"),
    ("Loop", "10s–100s kb", "held by CTCF + cohesin"),
    ("TAD", "0.1–1 Mb", "a neighbourhood of loops"),
    ("A/B compartment", "several Mb", "active with active"),
    ("Chromosome\nterritory", "a whole chromosome", "its own zone in the nucleus"),
]

def cartoon(k):
    ax = fig.add_axes([X[k], CY, W, CH]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal")
    ax.axis("off"); return ax

# 1 DNA
a = cartoon(0); t = np.linspace(0.08, 0.92, 200)
for ph, c in ((0, CRIMSON), (np.pi, DARK)):
    a.plot(t, 0.5 + 0.13 * np.sin(2 * np.pi * 3 * (t - 0.08) + ph), color=c, lw=1.6)
for x in np.linspace(0.12, 0.88, 14):
    y1 = 0.5 + 0.13 * np.sin(2 * np.pi * 3 * (x - 0.08)); y2 = 0.5 + 0.13 * np.sin(2 * np.pi * 3 * (x - 0.08) + np.pi)
    a.plot([x, x], [y1, y2], color=LINE, lw=0.8, zorder=0)
# 2 nucleosomes
a = cartoon(1); xs = np.linspace(0.05, 0.95, 200)
a.plot(xs, 0.5 + 0.06 * np.sin(xs * 22), color=TEXT2, lw=1.2, zorder=1)
for x in (0.2, 0.5, 0.8):
    a.add_patch(mp.Circle((x, 0.5), 0.1, fc=PROT, ec=TEXT2, lw=1, zorder=2))
    a.plot([x + 0.03, x + 0.1], [0.59, 0.72], color=CRIMSON, lw=1.2, zorder=3)          # histone tail
    a.add_patch(mp.Circle((x + 0.1, 0.74), 0.025, fc=CRIMSON, ec="none", zorder=3))      # mark
a.text(0.5, 0.2, "tails carry the marks", ha="center", va="center", fontsize=7.5, color=CRIMSON)
# 3 loop
a = cartoon(2); cx, cy, R = 0.5, 0.58, 0.25
th = np.linspace(np.deg2rad(250), np.deg2rad(-70), 200)
a.plot(cx + R * np.cos(th), cy + R * np.sin(th), color=TEXT2, lw=1.6)
for ang, nx in ((250, 0.46), (290, 0.54)):
    px, py = cx + R * np.cos(np.deg2rad(ang)), cy + R * np.sin(np.deg2rad(ang))
    a.plot([px, nx], [py, 0.2], color=TEXT2, lw=1.6)
a.plot([0.05, 0.46], [0.2, 0.2], color=TEXT2, lw=1.6); a.plot([0.54, 0.95], [0.2, 0.2], color=TEXT2, lw=1.6)
a.add_patch(mp.Ellipse((0.5, 0.31), 0.22, 0.1, fc="none", ec=DARK, lw=1.8))
a.add_patch(mp.Circle(((cx + R * np.cos(np.deg2rad(200))), cy + R * np.sin(np.deg2rad(200))), 0.035,
                      fc=CRIMSON, ec="none", zorder=3))
a.add_patch(mp.Circle(((cx + R * np.cos(np.deg2rad(-20))), cy + R * np.sin(np.deg2rad(-20))), 0.035,
                      fc=tint(0.5), ec="none", zorder=3))
# 4 TAD: a rosette of loops inside one domain
a = cartoon(3)
a.add_patch(mp.Circle((0.5, 0.52), 0.36, fc=tint(0.08), ec=LINE, lw=1, ls=(0, (3, 2))))
for k, ang in enumerate((90, 210, 330)):
    r = np.deg2rad(ang); ccx, ccy = 0.5 + 0.17 * np.cos(r), 0.52 + 0.17 * np.sin(r)
    tt = np.linspace(0, 2 * np.pi, 80)
    a.plot(ccx + 0.13 * np.cos(tt), ccy + 0.13 * np.sin(tt), color=TEXT2, lw=1.3)
a.add_patch(mp.Ellipse((0.5, 0.52), 0.09, 0.05, fc="none", ec=DARK, lw=1.5))
a.plot([0.05, 0.14], [0.52, 0.52], color=TEXT2, lw=1.3); a.plot([0.86, 0.95], [0.52, 0.52], color=TEXT2, lw=1.3)
# 5 compartments: A segments gather with A, B with B
a = cartoon(4)
a.add_patch(mp.Circle((0.27, 0.6), 0.2, fc=tint(0.15), ec="none"))
a.add_patch(mp.Circle((0.73, 0.42), 0.2, fc="#EEF0F3", ec="none"))
a.text(0.27, 0.86, "A", ha="center", va="center", fontsize=10, fontweight="bold", color=CRIMSON)
a.text(0.73, 0.14, "B", ha="center", va="center", fontsize=10, fontweight="bold", color=DARK)
pts = [(0.2, 0.62), (0.7, 0.45), (0.32, 0.55), (0.78, 0.38), (0.25, 0.7), (0.68, 0.34), (0.34, 0.66)]
for (x0, y0), (x1, y1) in zip(pts[:-1], pts[1:]):
    a.plot([x0, x1], [y0, y1], color=LINE, lw=1, zorder=1)
for i, (x, y) in enumerate(pts):
    a.add_patch(mp.Circle((x, y), 0.045, fc=CRIMSON if i % 2 == 0 else GREYF, ec="none", zorder=2))
# 6 territories
a = cartoon(5)
a.add_patch(mp.Circle((0.5, 0.5), 0.44, fc="white", ec=TEXT2, lw=1.2))
for (x, y, w, h, ang, c) in ((0.33, 0.62, 0.32, 0.22, 30, tint(0.55)), (0.66, 0.65, 0.26, 0.2, -20, "#D5D9DF"),
                             (0.38, 0.33, 0.28, 0.2, -35, "#C9CDD3"), (0.68, 0.36, 0.24, 0.22, 15, "#E3E6EB")):
    a.add_patch(mp.Ellipse((x, y), w, h, angle=ang, fc=c, ec="none"))
a.text(0.5, 0.14, "nucleus", ha="center", va="center", fontsize=7.5, color=TEXT2)

# labels under each cartoon
for k, (name, size, what) in enumerate(levels):
    t = fig.add_axes([X[k], 0.31, W + GAP * 0.8, 0.2]); t.axis("off"); t.set_xlim(0, 1); t.set_ylim(0, 1)
    t.text(0, 0.78, name, ha="left", va="bottom", fontsize=10.5, fontweight="bold", color=TEXT, linespacing=1.1)
    t.text(0, 0.7, size, ha="left", va="top", fontsize=9, fontweight="semibold", color=CRIMSON)
    t.text(0, 0.5, what, ha="left", va="top", fontsize=8.3, color=TEXT2)

# arrow: getting bigger
ar = fig.add_axes([0, 0.79, 1, 0.04]); ar.axis("off"); ar.set_xlim(0, 1); ar.set_ylim(0, 1)
ar.annotate("", (0.995, 0.5), (0.005, 0.5), arrowprops=dict(arrowstyle="-|>", color=LINE, lw=1.2, mutation_scale=12))
ar.text(0.5, 0.95, "each level folds the one before it  →  1 m of DNA fits in a 5 µm nucleus", ha="center",
        va="bottom", fontsize=9, color=TEXT2)

# Hi-C row
n = 40; i = np.arange(n); D = np.abs(i[:, None] - i[None, :]); base = (D + 1.0) ** -0.85
def mat(k, M, note=None):
    m = fig.add_axes([X[k], MY, W, MH])
    M = rng.poisson(M * 300) / 300 + 1e-4; L = np.log10(M)
    m.imshow(L, cmap=CM, interpolation="nearest", vmin=np.percentile(L, 30), vmax=np.percentile(L, 99.7))
    m.set_xticks([]); m.set_yticks([])
    for s in m.spines.values(): s.set_color(LINE)
    return m
for k in (0, 1):
    t = fig.add_axes([X[k], MY, W, MH]); t.axis("off"); t.set_xlim(0, 1); t.set_ylim(0, 1)
    t.add_patch(mp.Rectangle((0, 0), 1, 1, fc=FAINT, ec="none"))
    t.text(0.5, 0.5, "too small for Hi-C\n" + ("" if k == 0 else "(Micro-C reaches\nthis level)"),
           ha="center", va="center", fontsize=8.3, color=TEXT2, linespacing=1.3)
M = base.copy(); M[8, 30] = M[30, 8] = 0.5
m = mat(2, M); m.add_patch(mp.Circle((30, 8), 2.6, fc="none", ec=TEXT, lw=1))
M = base.copy()
for a0, b0 in ((0, 13), (13, 27), (27, 40)): M[a0:b0, a0:b0] *= 2.6
mat(3, M)
c = np.where((i // 7) % 2 == 0, 1, -1); M = base * (1 + 0.8 * np.outer(c, c)) + 0.02 * (np.outer(c, c) > 0)
mat(4, M)
M = np.full((n, n), 0.012); M[:20, :20] = base[:20, :20] + 0.05; M[20:, 20:] = base[:20, :20] + 0.05
m = mat(5, M)
m.text(-1.2, 10, "chr1", ha="right", va="center", fontsize=7.5, color=TEXT2, rotation=90)
m.text(-1.2, 30, "chr2", ha="right", va="center", fontsize=7.5, color=TEXT2, rotation=90)
cap = fig.add_axes([0, MY - 0.055, 1, 0.045]); cap.axis("off"); cap.set_xlim(0, 1); cap.set_ylim(0, 1)
for k, s in ((2, "dot off the diagonal"), (3, "blocks on the diagonal"), (4, "checkerboard"),
             (5, "cis bright, trans dim")):
    cap.text(X[k], 0.5, s, ha="left", va="center", fontsize=8.5, color=TEXT, fontweight="semibold")
lab = fig.add_axes([0, MY + MH + 0.012, 1, 0.03]); lab.axis("off"); lab.set_xlim(0, 1); lab.set_ylim(0, 1)
lab.text(0, 0.5, "In a Hi-C contact matrix", ha="left", va="center", fontsize=9.5, fontweight="semibold",
         color=TEXT)
lab.plot([0.17, 1], [0.5, 0.5], color="#E3E6EB", lw=0.9)

hd = fig.add_axes([0, 0.9, 1, 0.01]); hd.axis("off")
header(hd, "3D genome · folding the genome",
       "Nucleosomes are only the first fold; Hi-C sees the folds above them", dy=0)
fig.savefig(OUT + "genome_folding.png"); print("ok g29 genome folding")
