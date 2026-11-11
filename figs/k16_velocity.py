from gallery import *
rng = np.random.default_rng(16)

# RNA velocity: per gene, unspliced vs spliced. Above the steady-state line = being switched on; below = off.
fig = plt.figure(figsize=(9.8, 4.6))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.95], wspace=0.28, left=0.06, right=0.99, top=0.8, bottom=0.12)
DARK = "#5E6670"; GREYF = "#B8BEC7"
# simulate one gene: du/dt = a - b u, ds/dt = b u - g s; on for t < T, then off
a_on, b, g, T, dt = 2.0, 1.0, 0.4, 7.0, 0.01
u = s = 0.0; U, S, ON = [], [], []
for i in range(int(16 / dt)):
    a = a_on if i * dt < T else 0.0
    u, s = u + dt * (a - b * u), s + dt * (b * u - g * s); U.append(u); S.append(s); ON.append(a > 0)
U, S, ON = map(np.array, (U, S, ON))
k = rng.choice(len(U), 260, replace=False)
uu = U[k] + rng.normal(0, 0.06, 260); ss = S[k] + rng.normal(0, 0.12, 260)
p = fig.add_subplot(gs[0])
p.plot(S[ON], U[ON], color=CRIMSON, lw=1, alpha=0.6); p.plot(S[~ON], U[~ON], color=DARK, lw=1, alpha=0.6)
p.scatter(ss, uu, s=8, color=np.where(ON[k], tint(0.6), GREYF), lw=0)
xs = np.linspace(0, 5.3, 10); p.plot(xs, g / b * xs, color=TEXT, lw=1, ls=(0, (4, 3)))
p.set_xlim(-0.2, 5.6); p.set_ylim(-0.15, 2.4); fig.canvas.draw()
q1, q2 = p.transData.transform([(1, g / b), (4, 4 * g / b)])
p.text(4.5, g / b * 4.5 - 0.3, "steady state", ha="center", va="center", fontsize=8.2, color=TEXT,
       rotation=np.degrees(np.arctan2(q2[1] - q1[1], q2[0] - q1[0])), rotation_mode="anchor")
p.text(1.2, 1.3, "above the line:\nswitching ON\n(extra unspliced)", fontsize=8.3, color=CRIMSON,
       linespacing=1.15, va="center")
p.text(2.55, 0.72, "below:\nswitching OFF", fontsize=8.3, color=DARK, va="center", linespacing=1.15)
p.set_xlabel("spliced (mature) mRNA", color=TEXT2, fontsize=9); p.set_ylabel("unspliced (new) mRNA", color=TEXT2, fontsize=9)
p.set_xticks([]); p.set_yticks([])
for sp in ("left", "bottom"): p.spines[sp].set_color(LINE)
p.set_title("(a) One gene, every cell", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
p.text(0, 1.03, "the direction comes from this gap", transform=p.transAxes, fontsize=9, color=TEXT2)

# (b) arrows on an embedding
q = fig.add_subplot(gs[1])
t = np.sort(rng.uniform(0, 1, 70)); P = np.c_[-2 + 4 * t, 1.1 * np.sin(2.6 * t - 0.7)] + rng.normal(0, 0.13, (70, 2))
D = np.c_[np.full(70, 4.0), 1.1 * 2.6 * np.cos(2.6 * t - 0.7)]
D = D / np.linalg.norm(D, axis=1)[:, None] * 0.38 + rng.normal(0, 0.06, (70, 2))
q.scatter(P[:, 0], P[:, 1], s=9, color=GREYF, lw=0)
for (x, y), (dx, dy) in zip(P[::2], D[::2]):
    q.annotate("", (x + dx, y + dy), (x, y), arrowprops=dict(arrowstyle="-|>", color=CRIMSON, lw=0.8, mutation_scale=6))
q.set_xlim(-2.6, 2.8); q.set_ylim(-2.0, 2.0); q.set_xticks([]); q.set_yticks([])
for sp in q.spines.values(): sp.set_visible(True); sp.set_color("#E3E6EB")
q.set_title("(b) An arrow per cell", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
q.text(0, 1.03, "drawn on the UMAP, so it reads as motion", transform=q.transAxes, fontsize=9, color=TEXT2)

# (c) caveats
c = fig.add_subplot(gs[2]); c.axis("off"); c.set_xlim(0, 1); c.set_ylim(0, 1)
c.set_title("(c) Before believing the arrows", loc="left", fontsize=12.3, fontweight="semibold", color=CRIMSON, pad=22)
c.text(0, 1.03, "three assumptions that can fail", transform=c.transAxes, fontsize=9, color=TEXT2)
items = [("① constant rates", "splicing and decay rates are assumed\nfixed per gene; known to fail for some"),
         ("② the 2-D picture", "arrows inherit every UMAP distortion;\ndirection in 2-D ≠ direction in the data"),
         ("③ nuclei", "snRNA-seq is intron-rich by design;\nneeds a different model")]
for i, (h, d) in enumerate(items):
    y = 0.9 - i * 0.33
    c.text(0, y, h, ha="left", va="center", fontsize=10, fontweight="semibold", color=TEXT)
    c.text(0, y - 0.13, d, ha="left", va="center", fontsize=8.5, color=TEXT2, linespacing=1.25)
header(p, "Single-cell RNA-seq · RNA velocity", "Where the direction comes from, and why it is contested", dy=40, x1=3.55)
fig.savefig(OUT + "rna_velocity.png"); print("ok k16 velocity")
