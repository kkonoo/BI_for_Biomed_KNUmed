from gallery import *
rng = np.random.default_rng(13)

# CLR: (a) each count is compared with the sample's geometric mean, so sequencing depth cancels out;
# (b) zeros must be replaced by a pseudocount, and with many zeros the choice changes the answer (computed)
fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.3), gridspec_kw={"wspace": 0.32, "width_ratios": [1.25, 1]})
DARK = "#5E6670"; GREYF = "#B8BEC7"

a = axes[0]
x = np.array([3000, 800, 250, 20, 4.0]); names = ["A", "B", "C", "D", "E"]
gm = lambda v: np.exp(np.log(v).mean())
clr = np.log(x / gm(x))
for row, (scale, lab) in enumerate(((1, f"{int(x.sum()):,} reads"), (10, f"{int(10 * x.sum()):,} reads"))):
    y = 1 - row
    v = x * scale; g = gm(v)
    a.plot([2, 4e4], [y, y], color=LINE, lw=0.8, zorder=0)
    a.plot([g, g], [y - 0.2, y + 0.08], color=CRIMSON, lw=1.6)
    a.text(g, y - 0.26, f"geometric mean\n{g:.0f}", ha="center", va="top", fontsize=7.8, color=CRIMSON, linespacing=1.1)
    for k, (vv, nm, cc) in enumerate(zip(v, names, clr)):
        a.scatter([vv], [y], s=46, color=DARK, zorder=3, ec="white", lw=0.8)
        a.text(vv, y + 0.13, f"{nm}\n{cc:+.2f}", ha="center", va="bottom", fontsize=8.2, color=TEXT, linespacing=1.05)
    a.text(2, y + 0.42, lab, ha="left", va="center", fontsize=9, color=TEXT, fontweight="semibold")
a.set_xscale("log"); a.set_xlim(2, 4e4); a.set_ylim(-0.75, 1.6)
a.set_yticks([]); a.set_xlabel("read count (log scale)", color=TEXT2, fontsize=9)
for sp in ("left", "top", "right"): a.spines[sp].set_visible(False)
a.spines["bottom"].set_color(LINE); a.tick_params(colors=TEXT2, labelsize=8.3)
a.set_title("(a) Count ÷ geometric mean, then log", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
a.text(0, 1.03, "ten times the reads: every dot moves, the CLR numbers do not", transform=a.transAxes, fontsize=9,
       color=TEXT2)

# (b) pseudocount: one community sequenced shallow (2,000 reads) and deep (50,000 reads)
b = axes[1]
lp = rng.normal(0, 2.0, 200); p = np.exp(lp) / np.exp(lp).sum()
i = np.argsort(-p)[4]                                     # a common taxon
pcs = [1, 0.5, 0.1, 0.01]; fc = []; zs = []; zd = []
for c in pcs:
    d = []
    for r in range(50):
        deep = rng.multinomial(50000, p) + 0.0; sh = rng.multinomial(2000, p) + 0.0
        f = lambda v: np.log(v + c) - np.log(v + c).mean()
        d.append(f(sh)[i] - f(deep)[i]); zs.append((sh == 0).mean()); zd.append((deep == 0).mean())
    fc.append(np.exp(np.mean(d)))
Y = np.arange(len(pcs))[::-1]
for y, c, v in zip(Y, pcs, fc):
    b.barh(y, np.log2(v), height=0.55, color=CRIMSON if abs(np.log2(v)) > 0.3 else GREYF)
    b.text(np.log2(v) + (0.08 if v > 1 else -0.08), y, f"×{v:.2f}", ha="left" if v > 1 else "right", va="center",
           fontsize=9, fontweight="semibold", color=TEXT)
b.axvline(0, color=TEXT, lw=0.8)
b.set_yticks(Y); b.set_yticklabels([f"pseudocount {c:g}" for c in pcs], fontsize=8.8, color=TEXT)
b.set_xlim(-2.2, 2.2); b.set_xticks([-2, -1, 0, 1, 2]); b.set_xticklabels(["¼", "½", "1", "2", "4"])
b.tick_params(axis="y", length=0); b.tick_params(colors=TEXT2, labelsize=8.3)
for sp in ("top", "right", "left"): b.spines[sp].set_visible(False)
b.spines["bottom"].set_color(LINE)
b.set_xlabel("apparent change of one taxon (truth = ×1)", color=TEXT2, fontsize=9)
b.set_title("(b) Zeros need a pseudocount", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
b.text(0, 1.03, f"same community: {np.mean(zs):.0%} zeros at 2,000 reads, {np.mean(zd):.0%} at 50,000",
       transform=b.transAxes, fontsize=9, color=TEXT2)
fig.text(0.56, -0.05, "No pseudocount is right in general; 0.1 happens to cancel out in this example.", fontsize=8.6,
         color=TEXT2)
header(a, "Metagenomics · CLR", "CLR removes the sequencing depth — but zeros bring it back", dy=40, x1=2.3)
fig.savefig(OUT + "clr.png")
print("ok n13 clr", np.round(clr, 2), "gm", round(gm(x), 1), "fc", np.round(fc, 2), "zeros", round(np.mean(zs), 2),
      round(np.mean(zd), 2))
