from gallery import *
rng = np.random.default_rng(21)

# DDA / DIA / targeted: what gets fragmented at one moment (top), and which peptides end up measured
# in each of 10 runs (bottom, simulated)
fig = plt.figure(figsize=(9.8, 6.2))
gs = fig.add_gridspec(2, 3, height_ratios=[0.8, 1.25], hspace=0.62, wspace=0.16, left=0.02, right=0.98,
                      top=0.84, bottom=0.07)
DARK = "#5E6670"; GREYF = "#B8BEC7"

# --- simulation: 36 peptides in 6 co-elution groups of 6, 10 runs ---
P, R, G, TOPN = 36, 10, 6, 3
true = rng.normal(0, 1.1, P)
grp = rng.permutation(np.repeat(np.arange(P // G), G))
obs = true[:, None] + rng.normal(0, 0.6, (P, R))           # run-to-run intensity noise
dda = np.zeros((P, R), bool)
for r in range(R):
    for g in range(P // G):
        idx = np.where(grp == g)[0]
        pick = idx[np.argsort(-obs[idx, r])[:TOPN + rng.integers(0, 2)]]   # top-N, N varies slightly with cycle time
        dda[pick, r] = True
dda &= obs > -1.6
dia = obs > -1.6                                           # fragmented regardless; lost only below the limit
order = np.argsort(-true)
targets = order[[2, 9, 16, 23, 30]]                        # a predefined list of five peptides
tgt = np.zeros((P, R), bool); tgt[targets] = obs[targets] > -2.6                                  # targeted scans are the most sensitive

# --- top row: one MS1 snapshot ---
mz = np.linspace(400, 1000, 900)
peaks = [(438, 0.95), (471, 0.30), (505, 0.62), (552, 0.18), (588, 1.00), (631, 0.24), (668, 0.44),
         (702, 0.13), (745, 0.71), (788, 0.21), (824, 0.35), (869, 0.10), (902, 0.52), (946, 0.16)]
spec = sum(h * np.exp(-0.5 * ((mz - c) / 1.6) ** 2) for c, h in peaks)
modes = [("DDA", "the top-N most intense peaks"), ("DIA", "every window, in turn"), ("Targeted (PRM)", "only a predefined list")]
for j, (nm, sub) in enumerate(modes):
    a = fig.add_subplot(gs[0, j])
    a.fill_between(mz, 0, spec, color=GREYF, lw=0)
    if j == 0:
        for c, h in sorted(peaks, key=lambda t: -t[1])[:4]:
            a.plot([c], [h + 0.1], marker="v", color=CRIMSON, ms=6, mec="none")
    elif j == 1:
        for k, x0 in enumerate(np.arange(400, 1000, 50)):
            a.add_patch(mp.Rectangle((x0, 0), 50, 1.25, fc=tint(0.08 if k % 2 else 0.18), ec="none", zorder=0))
    else:
        for c in (471, 788):
            a.add_patch(mp.Rectangle((c - 6, 0), 12, 1.25, fc=tint(0.18), ec="none", zorder=0))
            a.plot([c], [1.13], marker="v", color=CRIMSON, ms=6, mec="none")
    a.set_xlim(400, 1000); a.set_ylim(0, 1.25); a.set_yticks([]); a.set_xticks([])
    for sp in ("left", "top", "right"): a.spines[sp].set_visible(False)
    a.spines["bottom"].set_color(LINE)
    a.set_xlabel("MS1  m/z", color=TEXT2, fontsize=8.5, labelpad=2)
    a.set_title(nm, loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=20)
    a.text(0, 1.04, "fragments " + sub, transform=a.transAxes, fontsize=9, color=TEXT2)
    if j == 0: top0 = a

# --- bottom row: measured / missing in 10 runs ---
mats = [(dda[order], "DDA"), (dia[order], "DIA"), (tgt[order], "Targeted")]
onlist = np.isin(order, targets)
for j, (M, nm) in enumerate(mats):
    a = fig.add_subplot(gs[1, j])
    for i in range(P):
        for r in range(R):
            if j == 2 and not onlist[i]:
                fc = FAINT
            else:
                fc = GREYF if M[i, r] else CRIMSON
            a.add_patch(mp.Rectangle((r + 0.06, P - 1 - i + 0.08), 0.88, 0.84, fc=fc, ec="none"))
    a.set_xlim(0, R); a.set_ylim(-0.2, P); a.axis("off")
    sel = M[onlist] if j == 2 else M
    miss = 100 * (1 - sel.mean())
    a.text(0, P + 0.8, f"{miss:.0f}% missing", ha="left", va="bottom", fontsize=10.5, fontweight="bold",
           color=CRIMSON if miss > 0 else TEXT)
    if j == 2:
        a.text(R, P + 0.8, f"but only {onlist.sum()} peptides", ha="right", va="bottom", fontsize=9, color=TEXT2)
    if j == 0:
        a.text(-0.2, P - 0.5, "high", ha="right", va="center", fontsize=8, color=TEXT2)
        a.text(-0.2, 0.5, "low", ha="right", va="center", fontsize=8, color=TEXT2)
    a.text(R / 2, -0.8, "10 runs →", ha="center", va="top", fontsize=8.5, color=TEXT2)
fig.text(0.02, 0.51, f"Which peptides are measured in 10 runs?  (simulated, {P} peptides sorted by abundance;"
         "  grey = measured, red = missing, pale = not on the list)", fontsize=9.5, color=TEXT, fontweight="semibold")
header(top0, "Proteomics · acquisition", "DDA picks by chance; DIA fragments everything", dy=34, x1=3.35)
fig.savefig(OUT + "dda_vs_dia.png")
print("ok m3", round(100 * (1 - dda.mean())), round(100 * (1 - dia.mean())), round(100 * (1 - tgt[targets].mean())))
