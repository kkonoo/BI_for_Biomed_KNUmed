from gallery import *
rng = np.random.default_rng(3)

# 16S vs shotgun: (a) what each approach reads from the same sample; (b) how far down each can name things
W, H = 9.8, 5.0
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, (H - 0.8) / H]); ax.set_xlim(0, W); ax.set_ylim(0, H - 0.8); ax.axis("off")
DARK = "#5E6670"; GREYF = "#B8BEC7"

def ptitle(x, t, sub):
    ax.text(x, H - 1.05, t, ha="left", va="center", fontsize=12.3, fontweight="semibold", color=TEXT)
    ax.text(x, H - 1.35, sub, ha="left", va="center", fontsize=9, color=TEXT2)

# (a) genomes in one sample: three bacteria + host; the 16S gene is a small block in each bacterial genome
X0 = 0.15; ptitle(X0, "(a) What gets sequenced", "one sample: three bacteria and the host")
gx, gw = X0 + 0.05, 2.1
genomes = [("bacterium 1", 2.95, DARK), ("bacterium 2", 2.6, DARK), ("bacterium 3", 2.25, DARK), ("host", 1.8, GREYF)]
for nm, y, c in genomes:
    ax.add_patch(mp.Rectangle((gx, y - 0.07), gw, 0.14, fc=c if nm == "host" else "#D5D9DF", ec="none"))
    if nm != "host":
        ax.add_patch(mp.Rectangle((gx + 0.45, y - 0.07), 0.16, 0.14, fc=CRIMSON, ec="none"))
ax.text(gx + 0.53, 3.17, "16S gene", ha="center", va="bottom", fontsize=8.2, color=CRIMSON)
ax.text(gx + gw + 0.08, 2.6, "bacterial\ngenomes", ha="left", va="center", fontsize=8.2, color=TEXT2, linespacing=1.1)
ax.text(gx + gw + 0.08, 1.8, "human DNA", ha="left", va="center", fontsize=8.2, color=TEXT2)
ax.plot([gx + gw + 0.03] * 2, [2.2, 3.0], color=LINE, lw=0.8)

# 16S: PCR copies of the 16S block only
cx = 3.55
ax.text(cx, 3.22, "16S amplicon", ha="left", va="center", fontsize=10, fontweight="bold", color=TEXT)
for k in range(12):
    ax.add_patch(mp.Rectangle((cx + (k % 6) * 0.27, 2.86 - (k // 6) * 0.2), 0.2, 0.1, fc=CRIMSON, ec="none"))
ax.text(cx, 2.38, "PCR copies one gene\n→ who is there (roughly)", ha="left", va="top", fontsize=8.4, color=TEXT2,
        linespacing=1.15)
# shotgun: random fragments from every genome, host included
ax.text(cx, 1.6, "shotgun", ha="left", va="center", fontsize=10, fontweight="bold", color=TEXT)
cols = [("#D5D9DF", 7), (CRIMSON, 1), (GREYF, 4)]
pool = [c for c, n in cols for _ in range(n)]
rng.shuffle(pool)
for k, c in enumerate(pool):
    w = rng.uniform(0.14, 0.26)
    ax.add_patch(mp.Rectangle((cx + (k % 6) * 0.27, 1.24 - (k // 6) * 0.2), w, 0.1, fc=c, ec="none"))
ax.text(cx, 0.76, "random pieces of all DNA\n→ who is there + which genes;\nhuman reads are wasted", ha="left",
        va="top", fontsize=8.4, color=TEXT2, linespacing=1.15)
ax.annotate("", (cx - 0.12, 2.85), (gx + gw + 0.75, 2.6),
            arrowprops=dict(arrowstyle="-|>", color=GREYF, lw=1, mutation_scale=8, connectionstyle="arc3,rad=-0.2"))
ax.annotate("", (cx - 0.12, 1.2), (gx + gw + 0.75, 1.8),
            arrowprops=dict(arrowstyle="-|>", color=GREYF, lw=1, mutation_scale=8, connectionstyle="arc3,rad=0.2"))

# (b) resolution ladder
X1 = 6.1; ptitle(X1, "(b) How far down can it name?", "typical; depends on region, database and depth")
levels = ["phylum", "family", "genus", "species", "strain", "gene function"]
marks = {"16S": [2, 2, 2, 1, 0, 0], "shotgun": [2, 2, 2, 2, 1, 2]}
cxs = {"16S": X1 + 1.75, "shotgun": X1 + 2.85}
for nm, x in cxs.items():
    ax.text(x, 3.22, nm, ha="center", va="center", fontsize=9.5, fontweight="bold", color=TEXT)
for i, lv in enumerate(levels):
    y = 2.82 - i * 0.43
    if lv == "gene function":
        ax.plot([X1, X1 + 3.4], [y + 0.215] * 2, color=LINE, lw=0.8)
    ax.text(X1 + 0.05, y, lv, ha="left", va="center", fontsize=9, color=TEXT)
    for nm, x in cxs.items():
        m = marks[nm][i]
        if m == 2:
            ax.add_patch(mp.Circle((x, y), 0.11, fc=DARK, ec="none"))
        elif m == 1:
            ax.add_patch(mp.Circle((x, y), 0.11, fc="white", ec=CRIMSON, lw=1.4))
            ax.add_patch(mp.Wedge((x, y), 0.11, 90, 270, fc=CRIMSON, ec="none"))
        else:
            ax.add_patch(mp.Circle((x, y), 0.11, fc="white", ec=LINE, lw=1.2))
ly = 0.15
for k, (lab, m) in enumerate((("reliable", 2), ("unreliable", 1), ("no", 0))):
    x = X1 + 0.05 + k * 1.15
    if m == 2: ax.add_patch(mp.Circle((x + 0.08, ly), 0.07, fc=DARK, ec="none"))
    elif m == 1:
        ax.add_patch(mp.Circle((x + 0.08, ly), 0.07, fc="white", ec=CRIMSON, lw=1.2))
        ax.add_patch(mp.Wedge((x + 0.08, ly), 0.07, 90, 270, fc=CRIMSON, ec="none"))
    else: ax.add_patch(mp.Circle((x + 0.08, ly), 0.07, fc="white", ec=LINE, lw=1.1))
    ax.text(x + 0.22, ly, lab, ha="left", va="center", fontsize=8.4, color=TEXT2)
header(ax, "Metagenomics · 16S vs shotgun", "16S reads one gene; shotgun reads every gene in the tube", x=0.15 / W,
       dy=0, x1=1 - 0.15 / W)
fig.savefig(OUT + "amplicon_vs_shotgun.png"); print("ok n3")
