from gallery import *
rng = np.random.default_rng(19)

# GSEA uses the whole ranked list (left); ORA cuts a list at a threshold and runs Fisher's test (right)
N = 2000
fig = plt.figure(figsize=(9.8, 5.4))
gs = fig.add_gridspec(3, 2, height_ratios=[2.4, 0.38, 1.3], width_ratios=[1.15, 1], hspace=0.12, wspace=0.16,
                      left=0.06, right=0.99, top=0.84, bottom=0.09)
DARK = "#5E6670"; GREYF = "#B8BEC7"

# (a) GSEA ----------------------------------------------------------------------------------------------
rank = np.arange(N)
inset = np.sort(np.concatenate([rng.choice(np.arange(0, 320), 34, replace=False),
                                rng.choice(np.arange(320, N), 16, replace=False)]))
hit = np.zeros(N, bool); hit[inset] = True
w = np.linspace(2.6, -2.6, N)
Nh = np.abs(w[hit]).sum()
run = np.cumsum(np.where(hit, np.abs(w) / Nh, -1 / (N - hit.sum())))
ax = fig.add_subplot(gs[0, 0])
ax.plot(rank, run, color=CRIMSON, lw=2)
ax.axhline(0, color=LINE, lw=1)
k = np.argmax(run)
ax.plot([k, k], [0, run[k]], color=TEXT, lw=1.2, ls=(0, (3, 2)))
ax.text(k + 30, run[k] + 0.03, f"enrichment score = {run[k]:.2f}", ha="left", va="bottom", fontsize=9.3, color=TEXT,
        fontweight="semibold")
ax.text(N * 0.98, 0.78, "walk down the list:\nstep up at each gene in the set,\nstep down otherwise",
        ha="right", va="top", fontsize=8.8, color=TEXT2, linespacing=1.4)
ax.set_ylim(-0.03, 0.82)
ax.set_ylabel("running score", color=TEXT2, fontsize=9)
ax.set_xlim(0, N); ax.set_xticks([])
for sp in ("left", "bottom"): ax.spines[sp].set_color(LINE)
ax.tick_params(colors=TEXT2, labelsize=8.5)
ax.set_title("(a) GSEA · the whole ranked list", loc="left", fontsize=12.5, fontweight="semibold", color=TEXT, pad=24)
ax.annotate("no threshold: every gene contributes", (0, 1), xycoords="axes fraction", xytext=(0, 6),
            textcoords="offset points", ha="left", va="bottom", fontsize=9.3, color=TEXT2)
axh = fig.add_subplot(gs[1, 0])
axh.vlines(inset, 0, 1, color=TEXT, lw=0.9)
axh.set_xlim(0, N); axh.set_ylim(0, 1); axh.axis("off")
axh.text(-15, 0.5, "genes in\nthe set", ha="right", va="center", fontsize=8.3, color=TEXT2, linespacing=1.2)
axr = fig.add_subplot(gs[2, 0])
axr.fill_between(rank, 0, w, where=w > 0, color=tint(0.45), lw=0)
axr.fill_between(rank, 0, w, where=w <= 0, color=GREYF, lw=0)
axr.axhline(0, color=LINE, lw=1)
axr.set_xlim(0, N); axr.set_yticks([]); axr.set_xticks([])
for sp in ("left", "bottom"): axr.spines[sp].set_color(LINE)
axr.text(30, -0.35, "up in treatment", ha="left", va="top", fontsize=8.5, color=CRIMSON)
axr.text(N - 30, 0.35, "down in treatment", ha="right", va="bottom", fontsize=8.5, color=DARK)
axr.set_xlabel("all 12,000 genes, ranked by the DE statistic", color=TEXT2, fontsize=9)

# (b) ORA -----------------------------------------------------------------------------------------------
b = fig.add_subplot(gs[:, 1]); b.set_xlim(0, 10); b.set_ylim(0, 10); b.axis("off")
b.set_title("(b) ORA · a cut-off list", loc="left", fontsize=12.5, fontweight="semibold", color=TEXT, pad=24)
b.annotate("Fisher's exact test on a 2×2 table", (0, 1), xycoords="axes fraction", xytext=(0, 6),
           textcoords="offset points", ha="left", va="bottom", fontsize=9.3, color=TEXT2)
# ① threshold cut
b.text(0, 9.65, "① cut the list", ha="left", va="center", fontsize=10, fontweight="semibold", color=TEXT)
b.add_patch(mp.Rectangle((0, 8.55), 9.6, 0.5, fc=GREYF, ec="none"))
b.add_patch(mp.Rectangle((0, 8.55), 0.25, 0.5, fc=CRIMSON, ec="none"))
b.text(9.55, 8.8, "12,000 expressed genes", ha="right", va="center", fontsize=8.8, color="white",
       fontweight="semibold")
b.annotate("64 pass padj < 0.05, |log₂FC| > 1", (0.12, 8.55), (0.6, 7.95), fontsize=8.8, color=CRIMSON,
           fontweight="semibold", ha="left", va="center",
           arrowprops=dict(arrowstyle="-", color=CRIMSON, lw=1, connectionstyle="angle,angleA=0,angleB=90,rad=0"))
b.text(0.6, 7.45, "the other 11,936 are thrown away", ha="left", va="center", fontsize=8.5, color=TEXT2)
# ② 2×2 table
b.text(0, 6.65, "② count the overlap with one pathway", ha="left", va="center", fontsize=10,
       fontweight="semibold", color=TEXT)
X0, Y0, CWd, RH = 3.2, 3.55, 2.9, 1.0
cols = ["significant", "not significant"]; rows = ["in pathway", "not in pathway"]
vals = [["12", "188"], ["52", "11,748"]]
for j, c in enumerate(cols):
    b.text(X0 + CWd * (j + 0.5), Y0 + 2 * RH + 0.35, c, ha="center", va="center", fontsize=8.8, color=TEXT2)
for i, r in enumerate(rows):
    y = Y0 + RH * (1 - i)
    b.text(X0 - 0.15, y + RH / 2, r, ha="right", va="center", fontsize=8.8, color=TEXT2)
    for j in range(2):
        hl = (i == 0 and j == 0)
        b.add_patch(mp.Rectangle((X0 + CWd * j, y), CWd, RH, fc=CRIMSON if hl else FAINT, ec="white", lw=2))
        b.text(X0 + CWd * (j + 0.5), y + RH / 2, vals[i][j], ha="center", va="center", fontsize=10.5,
               fontweight="bold" if hl else "normal", color="white" if hl else TEXT)
b.text(X0, Y0 - 0.55, "expected by chance ≈ 1  (64 × 200 / 12,000)\nobserved 12  →  about 11× enriched", ha="left",
       va="center", fontsize=8.8, color=TEXT, linespacing=1.4)
# ③ pitfalls
b.text(0, 2.0, "③ where it goes wrong", ha="left", va="center", fontsize=10, fontweight="semibold", color=TEXT)
for k, (t, d) in enumerate((("background", "use expressed genes, not all ~20,000"),
                            ("threshold", "0.05 vs 0.10 can change the top pathway"),
                            ("small changes", "50 genes × 1.3 each: missed by ORA, found by GSEA"))):
    y = 1.4 - k * 0.55
    b.text(0.15, y, "⚠ " + t, ha="left", va="center", fontsize=8.8, color=CRIMSON, fontweight="semibold")
    b.text(2.75, y, d, ha="left", va="center", fontsize=8.8, color=TEXT2)

header(ax, "Interpretation · gene set enrichment",
       "GSEA uses the whole ranking; ORA keeps only the genes above a cut-off", dy=44, x1=2.03)
fig.savefig(OUT + "ora_vs_gsea.png"); print("ok h6")
