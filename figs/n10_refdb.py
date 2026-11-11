from gallery import *

# Reference databases: (a) a read is named only as far down as its closest relative in the database;
# (b) names change between database versions (old name → current name)
W, H = 9.8, 4.6
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, (H - 0.8) / H]); ax.set_xlim(0, W); ax.set_ylim(0, H - 0.8); ax.axis("off")
DARK = "#5E6670"; GREYF = "#B8BEC7"

def ptitle(x, t, sub):
    ax.text(x, H - 1.05, t, ha="left", va="center", fontsize=12.3, fontweight="semibold", color=TEXT)
    ax.text(x, H - 1.35, sub, ha="left", va="center", fontsize=9, color=TEXT2)

X0 = 0.15; ptitle(X0, "(a) Named only as far as the database reaches", "filled = genome in the database")
# tree: leaves top to bottom; x = depth
ys = [2.95, 2.55, 2.05, 1.65, 1.05, 0.65]
inDB = [True, True, True, False, False, False]
lx = 2.05
def hline(x0, x1, y): ax.plot([x0, x1], [y, y], color=TEXT2, lw=1.1, solid_capstyle="butt")
def vline(x, y0, y1): ax.plot([x, x], [y0, y1], color=TEXT2, lw=1.1, solid_capstyle="butt")
# genus 1 = leaves 0,1 ; genus 2 = leaves 2,3 ; family-level novel clade = 4,5
g1, g2, g3 = 1.55, 1.55, 1.35
for i in (0, 1): hline(g1, lx, ys[i])
vline(g1, ys[1], ys[0])
for i in (2, 3): hline(g2, lx, ys[i])
vline(g2, ys[3], ys[2])
for i in (4, 5): hline(g3, lx, ys[i])
vline(g3, ys[5], ys[4])
f1 = 1.0
m12 = (ys[0] + ys[1]) / 2; m34 = (ys[2] + ys[3]) / 2; m56 = (ys[4] + ys[5]) / 2
hline(f1, g1, m12); hline(f1, g2, m34); vline(f1, m34, m12)
r = 0.45
mf = (m12 + m34) / 2
hline(r, f1, mf); hline(r, g3, m56); vline(r, m56, mf)
for y, d in zip(ys, inDB):
    ax.add_patch(mp.Circle((lx + 0.1, y), 0.085, fc=DARK if d else "white", ec=DARK, lw=1.2, zorder=3))
# where reads from each organism end up
res = [(ys[1], "read from here → named to species", DARK, "semibold"),
       (ys[3], "no genome of its own → named only to genus", CRIMSON, "semibold"),
       (m56, "no relative at all → 'unclassified'\n(microbial dark matter)", CRIMSON, "semibold")]
for y, t, c, wt in res:
    ax.text(lx + 0.32, y, t, ha="left", va="center", fontsize=8.7, color=c, fontweight=wt, linespacing=1.15)
ax.text(X0, 0.18, "Reference genomes come mostly from cultured bacteria and Western cohorts.", ha="left", va="center",
        fontsize=8.6, color=TEXT2)

X1 = 6.0; ptitle(X1, "(b) Same organism, new name", "grey = older name, black = current name")
pairs = [("Firmicutes", "Bacillota"), ("Clostridium difficile", "Clostridioides difficile"),
         ("Lactobacillus plantarum", "Lactiplantibacillus plantarum"),
         ("Bacteroides vulgatus", "Phocaeicola vulgatus"), ("Eubacterium rectale", "Agathobacter rectalis")]
for k, (o, n) in enumerate(pairs):
    y = 2.9 - k * 0.5
    ax.text(X1, y + 0.1, o, ha="left", va="center", fontsize=8.9, color=TEXT2)
    ax.text(X1 + 0.15, y - 0.13, "→ " + n, ha="left", va="center", fontsize=9.2, color=TEXT, fontweight="semibold")
ax.text(X1, 0.18, "Check which database and which version was used.", ha="left", va="center", fontsize=8.6,
        color=CRIMSON, fontweight="semibold")
header(ax, "Metagenomics · reference databases", "You can only name what the database already contains",
       x=0.15 / W, dy=0, x1=1 - 0.15 / W)
fig.savefig(OUT + "reference_db.png"); print("ok n10")
