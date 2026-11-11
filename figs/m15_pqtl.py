from gallery import *
rng = np.random.default_rng(15)

# pQTLs: (a) a simulated genome-wide cis/trans map; (b) the chain from GWAS hit to drug target
fig = plt.figure(figsize=(9.8, 5.0))
gs = fig.add_gridspec(1, 2, width_ratios=[1, 0.95], wspace=0.22, left=0.08, right=0.98, top=0.8, bottom=0.11)
DARK = "#5E6670"; GREYF = "#B8BEC7"
CHR = np.array([248, 242, 198, 190, 181, 171, 159, 145, 138, 134, 135, 133, 114, 107, 102, 90, 83, 80, 59, 64, 47, 51])
start = np.r_[0, np.cumsum(CHR)[:-1]]; G = CHR.sum()

a = fig.add_subplot(gs[0])
for k in range(0, 22, 2):
    a.axvspan(start[k], start[k] + CHR[k], color=FAINT, lw=0, zorder=0)
    a.axhspan(start[k], start[k] + CHR[k], color=FAINT, lw=0, zorder=0, alpha=0.6)
ng = 900
gene = rng.uniform(0, G, ng)
cis = gene + rng.normal(0, 2.0, ng)
a.scatter(cis, gene, s=5, color=DARK, lw=0, alpha=0.8, label="cis: variant near the protein's own gene")
tv = rng.uniform(0, G, 160); tg = rng.uniform(0, G, 160)
abo = start[8] + 133
hv = abo + rng.normal(0, 0.5, 70); hg = rng.uniform(0, G, 70)
a.scatter(np.r_[tv, hv], np.r_[tg, hg], s=7, color=CRIMSON, lw=0, alpha=0.75, label="trans: variant elsewhere")
a.text(abo + 60, G * 0.93, "one locus,\nmany proteins\n(e.g. ABO)", fontsize=8.5, color=CRIMSON, ha="left",
       va="center", linespacing=1.15, bbox=dict(fc="white", ec="none", pad=1.5), zorder=4)
a.set_xlim(0, G); a.set_ylim(0, G)
mid = start + CHR / 2
a.set_xticks(mid[[0, 4, 9, 14, 21]]); a.set_xticklabels(["1", "5", "10", "15", "22"])
a.set_yticks(mid[[0, 4, 9, 14, 21]]); a.set_yticklabels(["1", "5", "10", "15", "22"])
a.tick_params(colors=TEXT2, labelsize=8.3, length=0)
for sp in a.spines.values(): sp.set_color(LINE)
a.set_xlabel("variant position (chromosome)", color=TEXT2, fontsize=9)
a.set_ylabel("position of the protein's gene", color=TEXT2, fontsize=9)
a.legend(frameon=True, framealpha=1, edgecolor="none", fontsize=8.2, loc="lower right", labelcolor=TEXT,
         markerscale=2, handletextpad=0.2)
a.set_title("(a) Where pQTLs sit", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
a.text(0, 1.03, "simulated; one dot = one variant–protein association", transform=a.transAxes, fontsize=9, color=TEXT2)

# (b) chain
b = fig.add_subplot(gs[1]); b.set_xlim(0, 10); b.set_ylim(0, 10); b.axis("off")
steps = [("GWAS hit", "a locus is linked to disease", "W2"), ("pQTL", "the same locus changes a protein's level", "W6"),
         ("colocalisation", "one shared causal variant?", ""), ("Mendelian randomisation", "does the protein cause the disease?", "W2"),
         ("drug target", "with genetic support", "")]
for k, (h, d, w) in enumerate(steps):
    y = 9.0 - k * 2.0
    last = k == len(steps) - 1
    b.add_patch(mp.FancyBboxPatch((0.2, y - 0.62), 9.4, 1.24, boxstyle="round,pad=0,rounding_size=0.18",
                                  fc=CRIMSON if last else "white", ec=CRIMSON if last else LINE, lw=0.9))
    b.text(0.6, y + 0.2, h, ha="left", va="center", fontsize=10.5, fontweight="bold", color="white" if last else TEXT)
    b.text(0.6, y - 0.27, d, ha="left", va="center", fontsize=8.8, color="white" if last else TEXT2)
    if w: b.text(9.3, y, w, ha="right", va="center", fontsize=8.5, color=TEXT2)
    if not last:
        b.annotate("", (1.2, y - 1.38), (1.2, y - 0.64), arrowprops=dict(arrowstyle="-|>", color=GREYF, lw=1,
                                                                          mutation_scale=8))
b.text(1.6, 9.0 - 3.0, "cis-pQTL from an affinity assay? rule out an epitope effect", ha="left", va="center",
       fontsize=8.2, color=CRIMSON)
b.set_title("(b) From GWAS hit to drug target", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
b.text(0, 1.03, "each step answers one question", transform=b.transAxes, fontsize=9, color=TEXT2)
header(a, "Proteomics · pQTL", "Genetics links a disease locus to a protein — and the protein to a drug", dy=40,
       x1=2.12)
fig.savefig(OUT + "pqtl.png"); print("ok m15 pqtl")
