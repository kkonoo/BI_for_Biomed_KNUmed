from gallery import *

# 3' vs 5' kits: the barcode is attached at one end, and reads pile up at that end
W, H = 9.8, 5.5
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, (H - 0.8) / H]); ax.set_xlim(0, W); ax.set_ylim(0, H - 0.8); ax.axis("off")
DARK = "#5E6670"; MONO = "DejaVu Sans Mono"
X0, X1 = 0.95, 5.55
segs = [("5′ UTR", 0.95, 1.3, "#E3E6EB"), ("V(D)J", 1.3, 2.15, CRIMSON), ("C (constant)", 2.15, 3.75, tint(0.3)),
        ("3′ UTR", 3.75, 4.95, "#E3E6EB")]

def transcript(y):
    for nm, a, b, c in segs:
        ax.add_patch(mp.Rectangle((a, y - 0.07), b - a, 0.14, fc=c, ec="none", zorder=2))
        ax.text((a + b) / 2, y + 0.2, nm, ha="center", va="bottom", fontsize=8,
                color=CRIMSON if nm == "V(D)J" else TEXT2, fontweight="semibold" if nm == "V(D)J" else "normal")
    ax.text(5.0, y, "AAAA", ha="left", va="center", fontsize=7.8, color=TEXT2, family=MONO)
    ax.add_patch(mp.Circle((X0, y), 0.06, fc=DARK, ec="none", zorder=3))
    ax.text(X0 - 0.12, y, "5′", ha="right", va="center", fontsize=9, color=TEXT)
    ax.text(5.5, y, "3′", ha="left", va="center", fontsize=9, color=TEXT)

def bead(x, y, lab):
    ax.add_patch(mp.Circle((x, y), 0.22, fc=DARK, ec="none", zorder=3))
    ax.text(x, y, "bead", ha="center", va="center", fontsize=7.2, color="white", zorder=4)
    ax.text(x, y - 0.36, lab, ha="center", va="top", fontsize=7.4, color=TEXT2, family=MONO)

def coverage(y, three):
    xs = np.linspace(X0, 4.95, 200)
    f = np.exp(-(4.95 - xs) / 0.55) if three else np.exp(-(xs - X0) / 0.75)
    ax.fill_between(xs, y, y + 0.62 * f, color=tint(0.22), lw=0, zorder=1)
    ax.plot(xs, y + 0.62 * f, color=tint(0.7), lw=1.1, zorder=2)
    ax.plot([X0, 4.95], [y, y], color=LINE, lw=0.8)
    ax.text(0.15, y + 0.2, "where the\nreads land", ha="left", va="center", fontsize=8, color=TEXT2,
            linespacing=1.15)
    ax.add_patch(mp.Rectangle((1.3, y), 0.85, 0.66, fc="none", ec=CRIMSON, lw=0.9, ls=(0, (2, 2)), zorder=3))

def side(y, title, chem, how, good, bad):
    x = 6.75
    ax.text(x, y + 0.55, title, ha="left", va="center", fontsize=11.5, fontweight="bold", color=TEXT)
    ax.text(x, y + 0.22, chem, ha="left", va="center", fontsize=8.8, color=TEXT2)
    ax.text(x, y - 0.25, how, ha="left", va="center", fontsize=8.8, color=TEXT, linespacing=1.3)
    ax.text(x, y - 0.8, "✓  " + good, ha="left", va="center", fontsize=9, color=CRIMSON, fontweight="semibold")
    if bad:
        ax.text(x, y - 1.08, "✕  " + bad, ha="left", va="center", fontsize=9, color=DARK)

# 3' kit
y = 3.95
transcript(y); bead(5.98, y, "TTTT·UMI·BC")
ax.annotate("", (X0 + 0.05, y - 0.32), (4.95, y - 0.32),
            arrowprops=dict(arrowstyle="-|>", color=TEXT2, lw=0.9, mutation_scale=8))
ax.text(2.95, y - 0.45, "reverse transcription starts on the bead, at the tail", ha="center", va="top",
        fontsize=8, color=TEXT2)
coverage(2.65, True)
ax.text(1.72, 2.65 + 0.3, "no reads", ha="center", va="center", fontsize=8.3, color=CRIMSON, fontweight="semibold")
side(3.65, "3′ kit — the default", "bead oligo: barcode · UMI · poly(dT)",
     "poly(dT) on the bead grabs the poly(A) tail\n→ the barcode joins the 3′ end",
     "gene expression; most atlases", "V(D)J is too far from the tail")

ax.plot([0.1, W - 0.1], [2.15, 2.15], color="#E3E6EB", lw=0.9)

# 5' kit
y = 1.45
transcript(y); bead(0.5 - 0.12, y + 0.0, ""); ax.text(0.38, y - 0.36, "BC·UMI·TSO", ha="center", va="top",
                                                    fontsize=7.4, color=TEXT2, family=MONO)
ax.text(5.0, y - 0.25, "TTTT", ha="left", va="center", fontsize=7.4, color=TEXT2, family=MONO)
ax.annotate("", (X0 + 0.05, y - 0.32), (4.95, y - 0.32),
            arrowprops=dict(arrowstyle="-|>", color=TEXT2, lw=0.9, mutation_scale=8))
ax.text(2.95, y - 0.45, "RT from a free poly(dT) primer, then switches onto the bead's TSO at the cap",
        ha="center", va="top", fontsize=8, color=TEXT2)
coverage(0.15, False)
ax.text(1.72, 0.15 + 0.47, "covered", ha="center", va="center", fontsize=8.3, color=CRIMSON, fontweight="semibold")
side(1.15, "5′ kit — immune repertoire", "bead oligo: barcode · UMI · TSO",
     "the enzyme reaches the cap and copies the TSO\n→ the barcode joins the 5′ end",
     "TCR/BCR + gene expression, same cells", None)
header(ax, "Single-cell RNA-seq · 3′ vs 5′", "The end that gets the barcode is the end that gets read", dy=0)
fig.savefig(OUT + "three_vs_five_prime.png"); print("ok k7 3' vs 5'")
