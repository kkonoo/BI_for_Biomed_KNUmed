from gallery import *
from matplotlib.colors import LinearSegmentedColormap

# The 3C family: same proximity ligation, different part of the contact matrix read out
CM = LinearSegmentedColormap.from_list("hic", [tint(0.12), tint(0.35), tint(0.7), CRIMSON])
FW, FH = 9.8, 6.8
fig = plt.figure(figsize=(FW, FH))
EMPTY = "#F4F5F7"

def mat(n, mask, loops=()):
    i = np.arange(n); d = np.abs(i[:, None] - i[None, :])
    M = (d + 1.0) ** -0.9
    for a, b in loops:
        for x, y in ((a, b), (b, a)): M[x, y] += 0.35
    rgba = CM(np.clip(M / 0.6, 0, 1))
    rgba[~mask] = matplotlib.colors.to_rgba(EMPTY)
    return rgba

methods = [
    ("3C", "one vs one", "PCR with two primers:\nare i and j in contact?", "one"),
    ("4C", "one vs all", "start from one viewpoint,\nread everything it touches", "row"),
    ("5C", "many vs many", "every pair inside\na chosen region", "block"),
    ("Capture Hi-C", "chosen baits vs all", "probes pull out junctions\nthat contain a bait (promoters)", "baits"),
    ("Hi-C", "all vs all", "biotin pull-down of\nevery ligation junction", "all"),
    ("Micro-C", "all vs all, finer", "MNase instead of a restriction\nenzyme → nucleosome resolution", "fine"),
    ("HiChIP / ChIA-PET", "protein-anchored", "add a ChIP step: only contacts\nat CTCF / H3K27ac sites", "chip"),
]
n = 24; VP = 7; BAITS = (4, 12, 19); BOUND = (3, 9, 16, 21)
W = 0.19; GAPX = 0.062; MH = W * FW / FH
rows = [(0, 4), (4, 7)]
Y = [0.5, 0.06]
for r, (s0, s1) in enumerate(rows):
    k = s1 - s0
    x_start = 0.0
    for c, (name, q, how, kind) in enumerate(methods[s0:s1]):
        x0 = x_start + c * (W + GAPX)
        ax = fig.add_axes([x0, Y[r] + 0.1, W, MH * 0.82])
        nn = 48 if kind == "fine" else n
        I, J = np.meshgrid(np.arange(nn), np.arange(nn), indexing="ij")
        if kind == "one":
            mask = ((I == VP) & (J == 17)) | ((I == 17) & (J == VP))
        elif kind == "row":
            mask = (I == VP) | (J == VP)
        elif kind == "block":
            mask = (I >= 8) & (I <= 18) & (J >= 8) & (J <= 18)
        elif kind == "baits":
            mask = np.isin(I, BAITS) | np.isin(J, BAITS)
        elif kind == "chip":
            mask = np.isin(I, BOUND) | np.isin(J, BOUND)
        else:
            mask = np.ones((nn, nn), bool)
        loops = [(3, 9), (9, 16), (16, 21)] if kind in ("all", "chip") else []
        if kind == "fine":
            loops = [(6, 18), (18, 32), (32, 42)]
        img = mat(nn, mask, loops)
        if kind == "chip":                       # anchor-anchor contacts stand out most
            for a, b in [(3, 9), (9, 16), (16, 21), (3, 16), (9, 21)]:
                for x, y in ((a, b), (b, a)): img[x, y] = matplotlib.colors.to_rgba(CRIMSON)
        if kind == "one":
            for x, y in ((VP, 17), (17, VP)): img[x, y] = matplotlib.colors.to_rgba(CRIMSON)
        ax.imshow(img, interpolation="nearest")
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values(): sp.set_color(LINE); sp.set_visible(True)
        if kind == "one":
            ax.text(-0.8, VP, "i", ha="right", va="center", fontsize=9, color=CRIMSON, fontweight="bold")
            ax.text(17, -0.8, "j", ha="center", va="bottom", fontsize=9, color=CRIMSON, fontweight="bold")
        if kind == "row":
            ax.text(-0.8, VP, "view-\npoint", ha="right", va="center", fontsize=7.5, color=CRIMSON,
                    linespacing=1.1)
        if kind == "baits":
            for b in BAITS: ax.text(-0.8, b, "▸", ha="right", va="center", fontsize=8, color=CRIMSON)
        if kind == "chip":
            for b in BOUND: ax.text(-0.8, b, "●", ha="right", va="center", fontsize=6, color=CRIMSON)
        t = fig.add_axes([x0, Y[r] - 0.035, W + 0.03, 0.13]); t.axis("off")
        t.text(0, 1, name, ha="left", va="top", fontsize=11.5, fontweight="bold", color=TEXT)
        t.text(0, 0.68, q, ha="left", va="top", fontsize=9.5, fontweight="semibold", color=CRIMSON)
        t.text(0, 0.4, how, ha="left", va="top", fontsize=8.6, color=TEXT2, linespacing=1.3)

# legend in the empty fourth slot of the second row
lx0 = 3 * (W + GAPX)
lg = fig.add_axes([lx0, Y[1] + 0.1, W + 0.05, MH * 0.82]); lg.axis("off"); lg.set_xlim(0, 1); lg.set_ylim(0, 1)
lg.add_patch(mp.Rectangle((0.0, 0.83), 0.09, 0.07, fc=tint(0.5), ec="none"))
lg.text(0.13, 0.865, "measured part of\nthe contact matrix", ha="left", va="center", fontsize=9, color=TEXT2,
        linespacing=1.25)
lg.add_patch(mp.Rectangle((0.0, 0.6), 0.09, 0.07, fc=EMPTY, ec=LINE, lw=0.6))
lg.text(0.13, 0.635, "not measured", ha="left", va="center", fontsize=9, color=TEXT2)
lg.text(0.0, 0.42, "each square = the contact\nmatrix of one region\n(rows, columns = genome bins)",
        ha="left", va="top", fontsize=9, color=TEXT2, linespacing=1.35)

hd = fig.add_axes([0.0, 0.875, 1, 0.01]); hd.axis("off")
header(hd, "3D genome · the 3C family",
       "Same proximity ligation, different question: which pairs to read out", dy=0)
fig.savefig(OUT + "hic_3c_family.png"); print("ok g25 3c family")
