from gallery import *
rng = np.random.default_rng(5)

# Across genes, mRNA and protein correlate at r ≈ 0.4. (a) the scatter; (b) what that buys you:
# knowing a gene's mRNA level narrows the plausible range of its protein level only a little.
fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.6), gridspec_kw={"wspace": 0.32, "width_ratios": [1, 1.05]})
DARK = "#5E6670"; GREYF = "#B8BEC7"
n = 2500
x = rng.normal(0, 1, n)
y = 0.4 * x + np.sqrt(1 - 0.4 ** 2) * rng.normal(0, 1, n)
r = np.corrcoef(x, y)[0, 1]
SM, SP = 0.75, 0.9                                     # sd across genes in log10 units: mRNA, protein
m, p = 1 + SM * x, 3 + SP * y

a = axes[0]
a.scatter(m, p, s=6, color=GREYF, lw=0, alpha=0.7)
b = np.polyfit(m, p, 1); xs = np.linspace(m.min(), m.max(), 10)
a.plot(xs, np.polyval(b, xs), color=DARK, lw=1.3)
a.text(0.03, 0.97, f"r = {r:.2f}\nr² = {r ** 2:.2f}", transform=a.transAxes, ha="left", va="top", fontsize=11,
       fontweight="bold", color=TEXT, linespacing=1.3)
a.set_xlabel("mRNA level (log₁₀)", color=TEXT2, fontsize=9.3)
a.set_ylabel("protein level (log₁₀)", color=TEXT2, fontsize=9.3)
for sp in ("left", "bottom"): a.spines[sp].set_color(LINE)
a.tick_params(colors=TEXT2, labelsize=8.3)
a.set_title("(a) One dot = one gene", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
a.text(0, 1.03, f"{n:,} simulated genes, correlation set to 0.4", transform=a.transAxes, fontsize=9, color=TEXT2)

# (b) 95% range of protein level: for a random gene vs given its mRNA level
c = axes[1]
res = p - np.polyval(b, m)
w_all = 2 * 1.96 * p.std(); w_given = 2 * 1.96 * res.std()
rows = [("any gene", w_all, GREYF), ("given its mRNA level", w_given, CRIMSON)]
for k, (lab, w, col) in enumerate(rows):
    yy = 1 - k
    c.plot([-w / 2, w / 2], [yy, yy], color=col, lw=9, solid_capstyle="butt")
    for t in np.arange(-w / 2, w / 2 + 1e-9, 1.0):             # one tick per order of magnitude
        c.plot([t, t], [yy - 0.09, yy + 0.09], color="white", lw=1.2)
    c.text(-w_all / 2, yy + 0.22, lab, ha="left", va="bottom", fontsize=10, fontweight="semibold", color=TEXT)
    c.text(w_all / 2, yy + 0.22, f"≈ {round(10 ** w, -2):,.0f}-fold", ha="right", va="bottom", fontsize=10.5,
           fontweight="bold", color=CRIMSON if k else TEXT2)
c.set_xlim(-w_all / 2 - 0.2, w_all / 2 + 0.2); c.set_ylim(-0.75, 1.75); c.axis("off")
c.text(-w_all / 2, -0.5, "white tick = one order of magnitude (10×)", ha="left", va="center", fontsize=8.6, color=TEXT2)
c.set_title("(b) How well can you guess protein from mRNA?", loc="left", fontsize=12.3, fontweight="semibold",
            color=TEXT, pad=22)
c.text(0, 1.03, "range that holds 95% of genes' protein levels", transform=c.transAxes, fontsize=9, color=TEXT2)
header(a, "Proteomics · why measure protein", "Knowing the mRNA still leaves a thousand-fold guess", dy=40, x1=2.38)
fig.savefig(OUT + "mrna_vs_protein.png"); print("ok m1", round(r, 2), round(10 ** w_all), round(10 ** w_given))
