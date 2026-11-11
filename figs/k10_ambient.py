from gallery import *
rng = np.random.default_rng(10)

# Ambient RNA: free mRNA from burst cells ("soup") enters every droplet, including empty ones
W, H = 9.8, 4.6
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, (H - 0.8) / H]); ax.set_xlim(0, W); ax.set_ylim(0, H - 0.8); ax.axis("off")
DARK = "#5E6670"; TC = "#8A929C"; HEP = "#5E6670"; SOUP = CRIMSON

def squiggle(x, y, s=0.12, c=SOUP, lw=1.1):
    t = np.linspace(0, 1, 20); ax.plot(x + t * s * 1.6, y + 0.03 * np.sin(t * 2 * np.pi * 1.5), color=c, lw=lw,
                                      solid_capstyle="round")
def ptitle(x, t, sub):
    ax.text(x, 3.6, t, ha="left", va="center", fontsize=11.3, fontweight="semibold", color=TEXT)
    ax.text(x, 3.33, sub, ha="left", va="center", fontsize=8.6, color=TEXT2)

# (a) suspension
ptitle(0.15, "(a) The cell suspension", "a burst cell leaks its mRNA")
ax.add_patch(mp.FancyBboxPatch((0.15, 0.7), 2.75, 2.4, boxstyle="round,pad=0,rounding_size=0.08", fc=FAINT, ec="none"))
for x, y in [(0.6, 2.6), (2.4, 2.65), (0.75, 1.15), (2.3, 1.2)]:
    ax.add_patch(mp.Circle((x, y), 0.24, fc=HEP, ec="none"))
for x, y in [(1.4, 2.75), (2.0, 1.9), (0.45, 1.85), (1.55, 1.05)]:
    ax.add_patch(mp.Circle((x, y), 0.12, fc=TC, ec="none"))
ax.add_patch(mp.Circle((1.35, 1.9), 0.26, fc="none", ec=SOUP, lw=1, ls=(0, (2, 2))))
ax.text(1.35, 1.9, "burst", ha="center", va="center", fontsize=7.4, color=SOUP)
for _ in range(40):
    x, y = rng.uniform(0.3, 2.6), rng.uniform(0.85, 2.95)
    if min(np.hypot(x - cx, y - cy) for cx, cy in [(0.6, 2.6), (2.4, 2.65), (0.75, 1.15), (2.3, 1.2), (1.35, 1.9), (1.4, 2.75),
                                                   (2.0, 1.9), (0.45, 1.85), (1.55, 1.05)]) > 0.36:
        squiggle(x, y)
ax.text(0.15, 0.45, "the soup comes mostly from the most\nabundant or most fragile cell type", ha="left", va="center",
        fontsize=8.3, color=TEXT2, linespacing=1.2)
ax.annotate("", (3.35, 1.9), (3.0, 1.9), arrowprops=dict(arrowstyle="-|>", color=LINE, lw=1.1, mutation_scale=9))

# (b) droplets
ptitle(3.45, "(b) Every droplet gets soup", "own RNA + a little soup")
drops = [("T cell + soup", TC, 0.12), ("hepatocyte + soup", HEP, 0.2), ("no cell: soup only", None, 0)]
for k, (lab, c, r) in enumerate(drops):
    y = 2.7 - k * 0.92
    ax.add_patch(mp.Circle((3.85, y), 0.36, fc="white", ec=LINE, lw=1))
    if c: ax.add_patch(mp.Circle((3.8, y + 0.04), r, fc=c, ec="none"))
    ax.add_patch(mp.Circle((4.02, y - 0.17), 0.06, fc=DARK, ec="none"))   # bead
    for dx, dy in [(-0.2, -0.2), (0.12, 0.27), (-0.15, 0.27)]:
        squiggle(3.85 + dx - 0.05, y + dy, s=0.07)
    ax.text(4.35, y, lab, ha="left", va="center", fontsize=8.8, color=TEXT)
ax.annotate("", (6.05, 1.9), (5.75, 1.9), arrowprops=dict(arrowstyle="-|>", color=LINE, lw=1.1, mutation_scale=9))

# (c) what you measure
ptitle(6.25, "(c) What you measure", "grey = own RNA · red = soup")
genes = ["CD3E", "IL7R", "ALB", "APOA1"]
own = [np.array([0.8, 0.6, 0, 0]), np.array([0, 0, 0.95, 0.75]), np.zeros(4)]
soup = np.array([0, 0, 0.22, 0.16])
for k in range(3):
    y = 2.42 - k * 0.92
    for i in range(4):
        x = 6.3 + i * 0.58
        ax.add_patch(mp.Rectangle((x, y), 0.42, own[k][i] * 0.6, fc=TC if k == 0 else HEP, ec="none"))
        ax.add_patch(mp.Rectangle((x, y + own[k][i] * 0.6), 0.42, soup[i] * 0.6, fc=SOUP, ec="none"))
        if k == 2: ax.text(x + 0.21, y - 0.13, genes[i], ha="center", va="center", fontsize=7.4, color=TEXT2)
    ax.plot([6.25, 8.55], [y, y], color=LINE, lw=0.8)
ax.text(8.65, 2.55, "a T cell\n“expressing”\nALB", ha="left", va="center", fontsize=8.4, color=CRIMSON,
        fontweight="semibold", linespacing=1.15)
ax.text(8.65, 1.6, "barely\nchanged", ha="left", va="center", fontsize=8.4, color=TEXT2, linespacing=1.15)
ax.text(8.65, 0.68, "= the soup\nprofile", ha="left", va="center", fontsize=8.4, color=TEXT2, linespacing=1.15)
ax.text(3.45, 0.2, "Correction: learn the soup from empty droplets, subtract it from every cell — SoupX, CellBender, DecontX",
        ha="left", va="center", fontsize=8.8, color=TEXT)
header(ax, "Single-cell RNA-seq · ambient RNA", "The soup gets into every droplet, including the empty ones", dy=0)
fig.savefig(OUT + "ambient_rna.png"); print("ok k10 ambient")
