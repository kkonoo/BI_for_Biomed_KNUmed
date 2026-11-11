from gallery import *
rng = np.random.default_rng(10)

# Target–decoy FDR: (a) decoy hits show what wrong matches score like; (b) 1% FDR for peptides is not 1% for proteins
fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.4), gridspec_kw={"wspace": 0.3, "width_ratios": [1.15, 1]})
DARK = "#5E6670"; GREYF = "#B8BEC7"

# (a) score distributions: 30,000 spectra, 55% have their true peptide in the database
N = 30000
correct = rng.random(N) < 0.55
t_sc = np.where(correct, rng.normal(4.2, 1.1, N), rng.gumbel(0, 0.6, N))   # best target match
d_sc = rng.gumbel(0, 0.6, N)                                                # best decoy match (always wrong)
is_dec = d_sc > t_sc                                                        # concatenated search: keep the top one
top = np.maximum(t_sc, d_sc)
target, decoy = top[~is_dec], top[is_dec]
thr_grid = np.linspace(0, 6, 601)
fdr = np.array([(decoy >= t).sum() / max((target >= t).sum(), 1) for t in thr_grid])
thr = thr_grid[np.argmax(fdr <= 0.01)]
a = axes[0]
bins = np.linspace(-2, 8, 70)
a.hist(target, bins, color=GREYF, label="target hits (true + wrong)")
a.hist(decoy, bins, histtype="step", color=CRIMSON, lw=1.6, label="decoy hits (all wrong)")
ym = a.get_ylim()[1]; a.set_ylim(0, ym * 1.5)
a.plot([thr, thr], [0, ym * 1.48], color=TEXT, lw=1, ls=(0, (3, 2)))
nacc = (target >= thr).sum(); ndec = (decoy >= thr).sum()
a.text(thr + 0.15, ym * 1.47, f"cut-off: accept {nacc:,}\ndecoys above it: {ndec}\n→ FDR ≈ {100 * ndec / nacc:.0f}%",
       ha="left", va="top", fontsize=8.8, color=TEXT, linespacing=1.3)
a.set_xlabel("match score", color=TEXT2, fontsize=9); a.set_yticks([])
for sp in ("left", "top", "right"): a.spines[sp].set_visible(False)
a.spines["bottom"].set_color(LINE); a.tick_params(colors=TEXT2, labelsize=8.3)
a.legend(frameon=False, fontsize=8.4, loc="upper left", labelcolor=TEXT, handlelength=1.2, bbox_to_anchor=(-0.02, 1.0))
a.set_title("(a) Decoys show what a wrong match scores", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT,
            pad=22)
a.text(0, 1.03, "each spectrum searched against real + reversed sequences", transform=a.transAxes, fontsize=9,
       color=TEXT2)

# (b) peptide → protein: false peptides scatter over the whole database
DB = 20000
present = rng.choice(DB, 3000, replace=False)
npep = rng.geometric(0.18, 3000)
true_prot = np.repeat(present, npep)
nfalse = int(round(len(true_prot) * 0.01 / 0.99))
allp = np.r_[true_prot, rng.integers(0, DB, nfalse)]
prots, cnt = np.unique(allp, return_counts=True)
wrong = ~np.isin(prots, present)
vals = [("peptides", 1.0), ("proteins", 100 * wrong.mean()), ("proteins with\none peptide", 100 * wrong[cnt == 1].mean()),
        ("proteins with\n≥ 2 peptides", 100 * wrong[cnt >= 2].mean())]
b = axes[1]
X = np.arange(len(vals))
cols = [DARK, CRIMSON, CRIMSON, GREYF]
b.bar(X, [v for _, v in vals], color=cols, width=0.6)
for x, (_, v) in zip(X, vals):
    b.text(x, v + 0.5, f"{v:.0f}%" if v >= 1 else f"{v:.1f}%", ha="center", va="bottom", fontsize=10,
           fontweight="bold", color=TEXT)
b.set_xticks(X); b.set_xticklabels([n for n, _ in vals], fontsize=8.6, color=TEXT)
b.tick_params(axis="x", length=0); b.set_yticks([])
for sp in ("left", "top", "right"): b.spines[sp].set_visible(False)
b.spines["bottom"].set_color(LINE)
b.set_ylim(0, max(v for _, v in vals) * 1.25)
b.set_title("(b) 1% of peptides ≠ 1% of proteins", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
b.text(0, 1.03, f"share wrong (simulated: {len(allp):,} peptides at 1% FDR)", transform=b.transAxes, fontsize=9,
       color=TEXT2)
b.text(0.0, 0.98, "true peptides pile up on the same\nproteins; wrong ones each land on\na different protein",
       transform=b.transAxes, ha="left", va="top", fontsize=8.5, color=TEXT2, linespacing=1.25)
header(a, "Proteomics · identification", "Error rates are estimated with decoys — and grow from peptide to protein",
       dy=40, x1=2.18)
fig.savefig(OUT + "target_decoy.png")
print("ok m10", round(thr, 2), nacc, ndec, [round(v, 1) for _, v in vals])
