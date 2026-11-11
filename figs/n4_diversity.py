from gallery import *
from n_sim import two_groups, bray_curtis, pcoa, permanova, rarefied_richness
rng = np.random.default_rng(29)

# Diversity, three standard figures: (a) richness and evenness can disagree; (b) observed richness depends on
# sequencing depth; (c) PCoA of Bray–Curtis distances with PERMANOVA. All computed on simulated communities.
fig, axes = plt.subplots(1, 3, figsize=(9.8, 4.1), gridspec_kw={"wspace": 0.42})
DARK = "#5E6670"; GREYF = "#B8BEC7"

def style(ax, t, sub):
    ax.set_title(t, loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
    ax.text(0, 1.03, sub, transform=ax.transAxes, fontsize=9, color=TEXT2)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"): ax.spines[sp].set_color(LINE)
    ax.tick_params(colors=TEXT2, labelsize=8.3)

# (a) X: 40 taxa, one dominant; Y: 12 taxa, even
a = axes[0]
X = np.r_[0.8, np.full(39, 0.2 / 39)]; Y = np.full(12, 1 / 12)
H = lambda p: -(p * np.log(p)).sum()
for yb, v, col in ((1.0, X, DARK), (0.0, Y, CRIMSON)):
    left = 0
    for k, f in enumerate(v):
        a.barh(yb, f, left=left, height=0.42, color=col if k % 2 == 0 else (GREYF if col == DARK else tint(0.55)),
               edgecolor="white", lw=0.3)
        left += f
a.text(0, 1.3, f"sample X:  richness {len(X)},  Shannon {H(X):.2f}", ha="left", va="bottom", fontsize=8.6,
       color=DARK, fontweight="semibold")
a.text(0, 0.3, f"sample Y:  richness {len(Y)},  Shannon {H(Y):.2f}", ha="left", va="bottom", fontsize=8.6,
       color=CRIMSON, fontweight="semibold")
a.text(0, -0.55, "one block = one taxon, width = its share", ha="left", va="center", fontsize=8.3, color=TEXT2)
a.set_xlim(0, 1); a.set_ylim(-0.7, 1.6); a.set_yticks([])
a.set_xticks([0, 0.5, 1]); a.set_xticklabels(["0", "50%", "100%"])
a.spines["left"].set_visible(False)
style(a, "(a) Alpha", "richer is not always more diverse")

# (b) rarefaction: same community at two depths, plus a poorer community
b = axes[1]
logp = rng.normal(0, 1.6, 300); p_rich = np.exp(logp) / np.exp(logp).sum()
logq = rng.normal(0, 1.6, 80); p_poor = np.exp(logq) / np.exp(logq).sum()
deep = rng.multinomial(40000, p_rich); shallow = rng.multinomial(6000, p_rich); poor = rng.multinomial(40000, p_poor)
dd = np.unique(np.r_[np.linspace(1, 40000, 80).astype(int), 6000])
for cnt, col, lab in ((deep, DARK, "sample P (40,000 reads)"), (poor, GREYF, "sample R: poorer community")):
    b.plot(dd, rarefied_richness(cnt, dd), color=col, lw=2)
ds = dd[dd <= 6000]
b.plot(ds, rarefied_richness(shallow, ds), color=CRIMSON, lw=2)
obs = {k: (v > 0).sum() for k, v in (("P", deep), ("Q", shallow), ("R", poor))}
b.scatter([40000, 6000, 40000], [obs["P"], obs["Q"], obs["R"]], s=22, color=[DARK, CRIMSON, GREYF], zorder=4)
b.text(39000, obs["P"] - 12, f"P: {obs['P']} taxa", ha="right", va="top", fontsize=8.5, color=DARK, fontweight="semibold")
b.text(6800, obs["Q"] - 4, f"Q: {obs['Q']}\nsame community,\nfewer reads", ha="left", va="top", fontsize=8.5,
       color=CRIMSON, fontweight="semibold", linespacing=1.15)
b.text(39000, obs["R"] + 10, f"R: {obs['R']} taxa", ha="right", va="bottom", fontsize=8.5, color=TEXT2,
       fontweight="semibold")
b.axvline(6000, color=TEXT2, lw=0.7, ls=(0, (3, 2)))
b.set_xlim(0, 41000); b.set_ylim(0, obs["P"] * 1.12)
b.set_xticks([0, 20000, 40000]); b.set_xticklabels(["0", "20k", "40k"])
b.set_xlabel("reads", color=TEXT2, fontsize=9); b.set_ylabel("taxa observed", color=TEXT2, fontsize=9)
style(b, "(b) Depth", "dashed = rarefy everyone to 6,000 reads")

# (c) beta: PCoA of Bray–Curtis, PERMANOVA
c = axes[2]
C, g = two_groups(60, rng, 1.0, n_eff=40)
D = bray_curtis(C); xy, expl = pcoa(D)
F, R2, p = permanova(D, g, rng)
c.scatter(*xy[g == 0].T, s=12, color=GREYF, lw=0, label="controls")
c.scatter(*xy[g == 1].T, s=12, color=CRIMSON, lw=0, alpha=0.85, label="cases")
c.set_xticks([]); c.set_yticks([])
c.set_xlabel(f"PCoA 1 ({expl[0]:.0%})", color=TEXT2, fontsize=9)
c.set_ylabel(f"PCoA 2 ({expl[1]:.0%})", color=TEXT2, fontsize=9)
ptxt = "p = 0.001" if p <= 0.001 else f"p = {p:.3f}"
c.text(0.98, 0.03, f"PERMANOVA\n{ptxt}, R² = {R2:.2f}", transform=c.transAxes, ha="right", va="bottom", fontsize=8.6,
       color=CRIMSON, fontweight="semibold", linespacing=1.2, bbox=dict(fc="white", ec="none", pad=1.5))
style(c, "(c) Beta", "Bray–Curtis; grey = 60 controls, red = 60 cases")
header(a, "Metagenomics · diversity", "Every diversity number depends on a choice — say which one you made",
       dy=40, x1=3.5)
fig.savefig(OUT + "diversity.png")
print("ok n4", {k: int(v) for k, v in obs.items()}, "H", round(H(X), 2), round(H(Y), 2), "R2", round(R2, 3), "p", p)
