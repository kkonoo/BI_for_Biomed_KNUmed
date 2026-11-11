from gallery import *

# Affinity platforms: Olink (two antibodies + DNA oligos) and SomaScan (aptamer); both read binding, not sequence,
# so a variant in the binding site lowers the signal without lowering the protein
W, H = 9.8, 4.5
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, (H - 0.8) / H]); ax.set_xlim(0, W); ax.set_ylim(0, H - 0.8); ax.axis("off")
ax.set_aspect("equal")
DARK = "#5E6670"; GREYF = "#B8BEC7"

def protein(x, y, r=0.36, c=GREYF):
    ax.add_patch(mp.Ellipse((x, y), 2 * r * 1.25, 2 * r, fc=c, ec="none", zorder=2))

def antibody(tip, ang, L=0.5, c=DARK, lw=2.4):
    """Y-shaped antibody: the two arms reach toward `tip`; ang = direction from tip to the stem (degrees)."""
    t = np.radians(ang); jx, jy = tip[0] + 0.26 * np.cos(t), tip[1] + 0.26 * np.sin(t)
    for d in (-32, 32):
        u = np.radians(ang + 180 + d)
        ax.plot([jx, jx + 0.3 * np.cos(u)], [jy, jy + 0.3 * np.sin(u)], color=c, lw=lw, solid_capstyle="round",
                zorder=3)
    sx, sy = jx + L * np.cos(t), jy + L * np.sin(t)
    ax.plot([jx, sx], [jy, sy], color=c, lw=lw, solid_capstyle="round", zorder=3)
    return sx, sy

def ptitle(x, t, sub):
    ax.text(x, H - 1.05, t, ha="left", va="center", fontsize=12.3, fontweight="semibold", color=TEXT)
    ax.text(x, H - 1.35, sub, ha="left", va="center", fontsize=9, color=TEXT2)

# (a) Olink PEA
X = 0.15; ptitle(X, "(a) Olink (PEA)", "two antibodies, each with a DNA tag")
cx, cy = 1.55, 1.3
protein(cx, cy)
ax.text(cx, cy, "target", ha="center", va="center", fontsize=8.3, color="white", fontweight="semibold", zorder=4)
s1 = antibody((cx - 0.3, cy + 0.3), 125)
s2 = antibody((cx + 0.3, cy + 0.3), 55)
top = max(s1[1], s2[1]) + 0.35
ax.plot([s1[0], cx - 0.02], [s1[1], top], color=CRIMSON, lw=1.6, zorder=3)
ax.plot([s2[0], cx + 0.02], [s2[1], top], color=CRIMSON, lw=1.6, zorder=3)
ax.plot([cx - 0.16, cx + 0.16], [top + 0.05, top + 0.05], color=CRIMSON, lw=1.6, zorder=3)
s1 = (s1[0], top - 0.33)
ax.text(cx, s1[1] + 0.68, "tags pair only if both bind", ha="center", va="center", fontsize=8.5, color=CRIMSON)
ax.text(X, 0.38, "→ copy the paired DNA, count it\n(qPCR or sequencing)", ha="left", va="center", fontsize=8.6,
        color=TEXT2, linespacing=1.2)

# (b) SomaScan
X = 3.35; ptitle(X, "(b) SomaScan", "one modified DNA aptamer")
cx = 4.6
protein(cx, cy)
tt = np.linspace(0, 1, 80)
ax.plot(cx - 0.05 + 0.16 * np.sin(tt * 4 * np.pi), cy + 0.36 + tt * 0.9, color=CRIMSON, lw=2, zorder=3)
ax.text(cx + 0.35, cy + 1.05, "aptamer folds\naround its target", ha="left", va="center", fontsize=8.5, color=CRIMSON,
        linespacing=1.15)
ax.text(X, 0.38, "→ read the bound aptamers\non a DNA array", ha="left", va="center", fontsize=8.6, color=TEXT2,
        linespacing=1.2)

# (c) binding, not identity
X = 6.45; ptitle(X, "(c) A variant in the binding site", "same amount of protein, different signal")
for k, (var, sig) in enumerate(((False, 1.0), (True, 0.3))):
    cx = 7.0 + k * 1.55
    protein(cx, cy)
    if var:
        ax.plot([cx - 0.18], [cy + 0.36], marker="*", ms=12, color=CRIMSON, mec="white", zorder=4)
        antibody((cx - 0.3, cy + 0.85), 110, c=GREYF)
        ax.text(cx + 0.05, cy + 1.25, "cannot\nbind", ha="left", va="center", fontsize=8.3, color=TEXT2, linespacing=1.1)
        ax.text(cx, cy - 0.65, "missense\nvariant", ha="center", va="top", fontsize=8.3, color=CRIMSON,
                linespacing=1.1)
    else:
        antibody((cx - 0.18, cy + 0.33), 110)
        ax.text(cx, cy - 0.65, "reference", ha="center", va="top", fontsize=8.3, color=TEXT2)
    bx = cx + 0.62
    ax.add_patch(mp.Rectangle((bx - 0.09, cy - 0.4), 0.18, 1.3 * sig, fc=DARK if not var else CRIMSON, ec="none"))
    ax.text(bx, cy - 0.5, "signal", ha="center", va="top", fontsize=7.5, color=TEXT2)
ax.text(X, 0.2, "looks like less protein → a false cis-pQTL", ha="left", va="center", fontsize=8.8,
        color=CRIMSON, fontweight="semibold")
header(ax, "Proteomics · affinity platforms", "Affinity assays measure binding, not identity", x=0.15 / W, dy=0,
       x1=1 - 0.15 / W)
fig.savefig(OUT + "affinity.png"); print("ok m14 affinity")
