from gallery import *
from matplotlib.colors import LinearSegmentedColormap
rng = np.random.default_rng(11)

# WGCNA: genes → co-expression modules → one summary per module → link to a trait
NS = 40                                           # samples
sizes = [34, 30, 26, 30]
MCOL = [CRIMSON, tint(0.45), "#5E6670", "#B8BEC7"]
MNAME = ["M1", "M2", "M3", "M4"]
trait = rng.normal(0, 1, NS)
drivers = [0.85 * trait + 0.53 * rng.normal(0, 1, NS)] + [rng.normal(0, 1, NS) for _ in range(3)]
X, lab = [], []
for m, (n, d) in enumerate(zip(sizes, drivers)):
    for _ in range(n):
        load = rng.uniform(0.55, 0.9) * rng.choice([1, 1, 1, -1])
        X.append(load * d + np.sqrt(1 - load ** 2) * rng.normal(0, 1, NS)); lab.append(m)
X = np.array(X); lab = np.array(lab)
C = np.abs(np.corrcoef(X))

def eigengene(m):
    Z = X[lab == m]; Z = (Z - Z.mean(1, keepdims=True)) / Z.std(1, keepdims=True)
    u, s, vt = np.linalg.svd(Z, full_matrices=False)
    e = vt[0] * np.sign(np.corrcoef(vt[0], Z.mean(0))[0, 1]); return e

fig = plt.figure(figsize=(9.8, 4.9))
gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.0, 1.0], wspace=0.32, left=0.03, right=0.99, top=0.76,
                      bottom=0.13)
def panel(ax, t, sub):
    ax.set_title(t, loc="left", fontsize=12.5, fontweight="semibold", color=TEXT, pad=22)
    ax.text(0, 1.02, sub, transform=ax.transAxes, ha="left", va="bottom", fontsize=9.3, color=TEXT2)
def style(ax):
    for sp in ("left", "bottom"): ax.spines[sp].set_color(LINE)
    ax.tick_params(colors=TEXT2, labelsize=8.5)

# ① correlation heatmap with module bars
a = fig.add_subplot(gs[0, 0])
CM = LinearSegmentedColormap.from_list("w", ["#FFFFFF", tint(0.3), CRIMSON, "#5E1220"])
a.imshow(C, cmap=CM, vmin=0, vmax=1, interpolation="nearest")
edges = np.cumsum([0] + sizes)
for m in range(4):
    a.add_patch(mp.Rectangle((edges[m] - 0.5, -6.5), sizes[m], 4.5, fc=MCOL[m], ec="white", lw=0.8,
                             clip_on=False))
    a.add_patch(mp.Rectangle((-6.5, edges[m] - 0.5), 4.5, sizes[m], fc=MCOL[m], ec="white", lw=0.8,
                             clip_on=False))
    a.text(-8, edges[m] + sizes[m] / 2 - 0.5, MNAME[m], ha="right", va="center", fontsize=8.5,
           color=TEXT, fontweight="semibold")
a.set_xticks([]); a.set_yticks([])
for sp in a.spines.values(): sp.set_visible(False)
a.set_xlim(-7, len(lab)); a.set_ylim(len(lab), -7)
a.set_xlabel("genes, ordered by module", color=TEXT2, fontsize=9)
panel(a, "① Genes that move together", "|correlation| between genes across samples")

# ② module eigengenes vs trait
b = fig.add_subplot(gs[0, 1]); style(b)
panel(b, "② One summary per module", "module eigengene vs a clinical trait")
for m, (c, mk) in enumerate(((CRIMSON, "o"), (MCOL[2], "o"))):
    mm = 0 if m == 0 else 2
    e = eigengene(mm); r = np.corrcoef(e, trait)[0, 1]
    b.scatter(trait, e, s=16, color=c, lw=0, alpha=0.85, zorder=3 - m)
    k, q = np.polyfit(trait, e, 1); xs = np.linspace(trait.min(), trait.max(), 10)
    b.plot(xs, k * xs + q, color=c, lw=1.4)
    b.text(0.03, 0.95 - 0.09 * m, f"{MNAME[mm]}: r = {r:.2f}", transform=b.transAxes, ha="left", va="top",
           fontsize=9.5, color=c, fontweight="semibold")
b.set_xlabel("trait (e.g. disease severity)", color=TEXT2, fontsize=9)
b.set_ylabel("module eigengene", color=TEXT2, fontsize=9)
b.set_xticks([]); b.set_yticks([])

# ③ why many samples: spread of r between two unrelated genes
c = fig.add_subplot(gs[0, 2]); style(c)
panel(c, "③ Needs many samples", "chance r between unrelated genes")
rr = np.linspace(-1, 1, 400)
for n, col, lw in ((6, "#B8BEC7", 1.6), (20, tint(0.55), 1.8), (50, CRIMSON, 2.0)):
    dens = (1 - rr ** 2) ** ((n - 4) / 2); dens /= dens.max()
    c.plot(rr, dens, color=col, lw=lw)
    xm = {6: 0.62, 20: 0.33, 50: 0.1}[n]
    c.text(xm + 0.05, (1 - xm ** 2) ** ((n - 4) / 2) + 0.03, f"n = {n}", ha="left", va="bottom", fontsize=9,
           color=col if n != 6 else "#5E6670", fontweight="semibold")
c.set_xlim(-1, 1); c.set_ylim(0, 1.15); c.set_yticks([])
c.set_xticks([-1, -0.5, 0, 0.5, 1]); c.set_xlabel("small n → big chance r → fake modules", color=TEXT2, fontsize=9)

a.apply_aspect(); pa = a.get_position()
for ax_ in (b, c):
    p_ = ax_.get_position(); ax_.set_position([p_.x0, pa.y0, p_.width, pa.height])
fig.text(0.03, 0.03, "20,000 genes → about a dozen modules → test each module summary against the trait  ·  "
         "⚠️ at least ~20 samples, better 50+", ha="left", va="center", fontsize=9.3, color=TEXT2)
header(a, "Beyond DE · co-expression network (WGCNA)",
       "Group genes that move together, then ask which group follows the trait", dy=46, x1=3.3)
fig.savefig(OUT + "wgcna.png"); print("ok h11 wgcna")
