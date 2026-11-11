from gallery import *

# Barcode and UMI: reads are grouped by cell barcode, then PCR copies sharing a UMI collapse to one molecule
W, H = 9.8, 4.3
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, (H - 0.8) / H]); ax.set_xlim(0, W); ax.set_ylim(0, H - 0.8); ax.axis("off")
DARK = "#5E6670"; MONO = "DejaVu Sans Mono"
# (a) what one read carries
ax.text(0.1, 3.3, "(a) What one read carries", ha="left", va="center", fontsize=11.6, fontweight="semibold", color=TEXT)
parts = [("cell barcode", "16 nt", DARK, 1.6), ("UMI", "12 nt", CRIMSON, 1.15), ("poly(dT)", "", "#C9CED6", 0.8),
         ("cDNA → which gene", "", tint(0.3), 2.4)]
x = 0.1
for nm, ln, c, w in parts:
    ax.add_patch(mp.Rectangle((x, 2.55), w, 0.36, fc=c, ec="white", lw=1.5))
    ax.text(x + w / 2, 2.73, nm, ha="center", va="center", fontsize=8.6, color="white" if c in (DARK, CRIMSON) else TEXT,
            fontweight="semibold")
    if ln: ax.text(x + w / 2, 2.35, ln, ha="center", va="center", fontsize=7.8, color=TEXT2)
    x += w
ax.text(0.1, 1.95, "cell barcode = which droplet (cell)", ha="left", va="center", fontsize=8.8, color=TEXT)
ax.text(0.1, 1.65, "UMI = which original molecule, tagged before PCR", ha="left", va="center", fontsize=8.8, color=CRIMSON)
ax.text(0.1, 1.35, "cDNA = which gene", ha="left", va="center", fontsize=8.8, color=TEXT)
ax.plot([6.35, 6.35], [0.2, 3.5], color="#E3E6EB", lw=0.9)
ax.plot([0.1, 6.1], [1.0, 1.0], color="#E3E6EB", lw=0.9)
# (b) collapsing reads
ax.text(0.1, 0.72, "(b) 7 reads, cell AAACCTG, gene X", ha="left", va="center", fontsize=10.5, fontweight="semibold",
        color=TEXT)
umis = [("GATTCAGG", CRIMSON, 3), ("TTAGCCAT", DARK, 2), ("CCGAATTC", tint(0.45), 2)]
x = 0.15
for u, c, n in umis:
    x1 = x
    for k in range(n):
        ax.add_patch(mp.Rectangle((x, 0.2), 0.42, 0.2, fc=c, ec="none")); x += 0.5
    ax.text((x1 + x - 0.08) / 2, 0.02, f"UMI {u[:5]}…", ha="center", va="center", fontsize=7.2, color=TEXT2, family=MONO)
    x += 0.2
ax.annotate("", (4.7, 0.3), (4.15, 0.3), arrowprops=dict(arrowstyle="-|>", color=TEXT2, lw=1, mutation_scale=9))
for i, (u, c, n) in enumerate(umis):
    ax.add_patch(mp.Circle((4.95 + i * 0.36, 0.3), 0.12, fc=c, ec="none"))
ax.text(5.95, 0.3, "= 3", ha="left", va="center", fontsize=11, fontweight="bold", color=TEXT)
ax.text(0.15, -0.32, "same UMI = PCR copies of one molecule → counted once", ha="left", va="center", fontsize=8.4, color=TEXT2)
# (c) the matrix
x0 = 6.65
ax.text(x0, 3.3, "(c) The count matrix", ha="left", va="center", fontsize=11.6, fontweight="semibold", color=TEXT)
ax.text(x0, 3.03, "molecules per cell × gene", ha="left", va="center", fontsize=8.7, color=TEXT2)
genes = ["gene X", "gene Y", "gene Z"]; cells = ["AAACCTG", "AAAGTAG", "AACTCCC", "…"]
vals = [[3, 0, 1], [0, 2, 0], [1, 0, 0], ["", "", ""]]
for j, g in enumerate(genes): ax.text(x0 + 1.35 + j * 0.62, 2.55, g, ha="center", va="center", fontsize=8, color=TEXT2)
for i, c in enumerate(cells):
    y = 2.15 - i * 0.42
    ax.text(x0, y, c, ha="left", va="center", fontsize=8, color=TEXT2, family=MONO)
    for j in range(3):
        v = vals[i][j]; hl = i == 0 and j == 0
        ax.add_patch(mp.Rectangle((x0 + 1.08 + j * 0.62, y - 0.17), 0.54, 0.34, fc=tint(0.15) if hl else FAINT,
                                  ec=CRIMSON if hl else "none", lw=1))
        ax.text(x0 + 1.35 + j * 0.62, y, str(v), ha="center", va="center", fontsize=9.5,
                color=CRIMSON if hl else (TEXT if v else "#9AA1AB"), fontweight="bold" if hl else "normal")
ax.text(x0, 0.3, "a count = molecules captured,\nnot reads sequenced", ha="left", va="center", fontsize=8.6,
        color=TEXT, linespacing=1.25)
header(ax, "Single-cell RNA-seq · barcodes and UMIs", "Every read says which cell and which molecule it came from", dy=0)
fig.savefig(OUT + "barcode_umi.png"); print("ok k20 umi")
