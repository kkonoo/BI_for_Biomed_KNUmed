from gallery import *
rng = np.random.default_rng(91)

# OTU vs ASV on one simulated amplicon sample. Three real sequences: A, B (2 bases from A) and C (14 bases from A),
# 253 bp (a V4 amplicon). Each read gets random substitution errors (0.4% per base). For every distinct sequence,
# plot how many times it was read against its distance from A.
Lseq, err = 253, 0.004
A = rng.integers(0, 4, Lseq)
def mutate(s, k):
    s = s.copy(); pos = rng.choice(Lseq, k, replace=False)
    s[pos] = (s[pos] + rng.integers(1, 4, k)) % 4; return s
B = mutate(A, 2); C = mutate(A, 14)
nreads = 3000
src = rng.choice(3, nreads, p=[0.68, 0.24, 0.08])
reads = np.array([(A, B, C)[k] for k in src])
E = rng.random(reads.shape) < err
reads[E] = (reads[E] + rng.integers(1, 4, E.sum())) % 4
uniq, counts = np.unique(reads, axis=0, return_counts=True)
dA = (uniq != A).sum(1)
is_true = [(u == A).all() or (u == B).all() or (u == C).all() for u in uniq]

# greedy 97% OTU clustering (abundance-sorted centroids, ≤ 7 mismatches = ≥ 97% identity)
order = np.argsort(-counts); cents = []
for i in order:
    if not any((uniq[i] != uniq[c]).sum() <= 7 for c in cents):
        cents.append(i)
n_otu = len(cents)

fig, (ax, side) = plt.subplots(1, 2, figsize=(9.8, 4.4), gridspec_kw={"width_ratios": [1.9, 1], "wspace": 0.12})
DARK = "#5E6670"; GREYF = "#B8BEC7"
ax.axvspan(-0.5, 7.5, color=FAINT, lw=0, zorder=0)
ax.text(5.0, 40, "one 97% OTU\n(≤ 7 differences\nfrom its centre)", ha="center", va="center", fontsize=8.6,
        color=TEXT2, linespacing=1.15)
x = dA + rng.uniform(-0.25, 0.25, len(dA))
y = counts * np.exp(rng.uniform(-0.12, 0.12, len(counts)))
tr = np.array(is_true)
ax.scatter(x[~tr], y[~tr], s=9, color=GREYF, lw=0, alpha=0.8, zorder=2)
for u, c_, lab in ((A, None, "A"), (B, None, "B"), (C, None, "C")):
    i = np.where((uniq == u).all(1))[0][0]
    ax.scatter([dA[i]], [counts[i]], s=46, color=CRIMSON, zorder=4, ec="white", lw=0.8)
    ax.text(dA[i] + 0.45, counts[i], f"{lab}  ×{counts[i]}", ha="left", va="center", fontsize=9.2,
            fontweight="semibold", color=CRIMSON)
ax.set_yscale("log"); ax.set_ylim(0.7, 2500); ax.set_xlim(-0.8, 17.5)
ax.set_xlabel("bases different from A (out of 253)", color=TEXT2, fontsize=9)
ax.set_ylabel("times this exact sequence was read", color=TEXT2, fontsize=9)
ax.set_xticks(range(0, 18, 2)); ax.set_yticks([1, 10, 100, 1000]); ax.set_yticklabels(["1", "10", "100", "1,000"])
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
for sp in ("left", "bottom"): ax.spines[sp].set_color(LINE)
ax.tick_params(colors=TEXT2, labelsize=8.6)
ax.text(11.3, 9, "grey = reads with\nsequencing errors:\neach version seen\nonly a few times", ha="left", va="center",
        fontsize=8.4, color=TEXT2, linespacing=1.15)
ax.set_title("One sample, three real sequences", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
ax.text(0, 1.03, f"simulated: {nreads:,} reads, 0.4% error per base, {len(uniq):,} distinct sequences",
        transform=ax.transAxes, fontsize=9, color=TEXT2)

side.axis("off"); side.set_xlim(0, 1); side.set_ylim(0, 1)
blocks = [("97% OTUs", f"{n_otu} units", "A and B fall inside one cluster\n→ merged; the cut-off was chosen to\nabsorb errors, "
           "and it absorbs B too", DARK),
          ("ASVs (DADA2, Deblur)", "3 units", "B is read too often, identically,\nto be an error of A\n→ kept as its own "
           "exact sequence", CRIMSON)]
for k, (h, n_, d, col) in enumerate(blocks):
    y0 = 0.93 - k * 0.47
    side.text(0, y0, h, ha="left", va="top", fontsize=10.3, fontweight="bold", color=TEXT)
    side.text(1, y0, n_, ha="right", va="top", fontsize=10.3, fontweight="bold", color=col)
    side.text(0, y0 - 0.1, d, ha="left", va="top", fontsize=8.6, color=TEXT2, linespacing=1.25)
side.plot([0, 1], [0.48, 0.48], color=LINE, lw=0.8)
header(ax, "Metagenomics · OTU vs ASV", "An error appears a few times; a real variant appears again and again",
       dy=40, x1=1.62)
fig.savefig(OUT + "otu_asv.png")
print("ok n9 uniq", len(uniq), "otus", n_otu, "A,B,C counts", [counts[np.where((uniq == u).all(1))[0][0]] for u in (A, B, C)],
      "max error count", counts[~tr].max())
