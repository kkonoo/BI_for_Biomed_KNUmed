from gallery import *
rng = np.random.default_rng(9)

# The human microbiome: (a) body sites are different ecosystems (typical phylum profiles, illustrative);
# (b) healthy guts differ a lot from person to person (simulated)
fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.6), gridspec_kw={"wspace": 0.22, "width_ratios": [1, 1.1]})
DARK = "#5E6670"; GREYF = "#B8BEC7"
phyla = [("Bacteroidetes", DARK), ("Firmicutes", CRIMSON), ("Proteobacteria", tint(0.45)),
         ("Actinobacteria", GREYF), ("Fusobacteria", "#8C939C"), ("other", "#E3E6EB")]
sites = {"Gut": [0.42, 0.45, 0.04, 0.04, 0.005, 0.045],
         "Mouth": [0.18, 0.36, 0.20, 0.12, 0.08, 0.06],
         "Skin": [0.05, 0.25, 0.15, 0.50, 0.005, 0.045],
         "Vagina": [0.02, 0.90, 0.02, 0.04, 0.0, 0.02],
         "Nose": [0.03, 0.32, 0.13, 0.48, 0.0, 0.04]}

def style(ax, t, sub):
    ax.set_title(t, loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
    ax.text(0, 1.03, sub, transform=ax.transAxes, fontsize=9, color=TEXT2)
    for sp in ("left", "top", "right"): ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_color(LINE); ax.tick_params(colors=TEXT2, labelsize=8.6, length=0)
    ax.set_yticks([0, 0.5, 1]); ax.set_yticklabels(["0", "50%", "100%"]); ax.set_ylim(0, 1)

def stack(ax, X, V, w, label_min=None):
    bl = np.zeros(len(X))
    for k, (nm, c) in enumerate(phyla):
        ax.bar(X, V[:, k], bottom=bl, color=c, width=w, edgecolor="white", lw=0.8)
        if label_min is not None:
            for x, v, b in zip(X, V[:, k], bl):
                if v >= label_min:
                    ax.text(x, b + v / 2, f"{v:.0%}", ha="center", va="center", fontsize=8,
                            color="white" if c in (DARK, CRIMSON, "#8C939C") else TEXT)
        bl += V[:, k]

a = axes[0]
V = np.array(list(sites.values()))
stack(a, np.arange(5), V, 0.62, label_min=0.1)
a.set_xticks(range(5)); a.set_xticklabels(sites.keys(), fontsize=9.3, color=TEXT)
style(a, "(a) Five body sites", "typical make-up by phylum (illustrative)")

b = axes[1]
n = 14
base = np.log(np.array(sites["Gut"]) + 0.005)
G = np.exp(base + rng.normal(0, 0.6, (n, 6))); G /= G.sum(1, keepdims=True)
G = G[np.argsort(-G[:, 0])]
stack(b, np.arange(n), G, 0.8)
b.set_xticks([]); b.set_xlabel("14 healthy adults (stool)", color=TEXT2, fontsize=9)
lo, hi = G[:, 0].min(), G[:, 0].max()
style(b, "(b) One body site, fourteen people", f"simulated; Bacteroidetes alone ranges from {lo:.0%} to {hi:.0%}")
handles = [mp.Patch(fc=c, ec="none", label=nm) for nm, c in phyla]
fig.legend(handles=handles, loc="lower center", ncol=6, frameon=False, fontsize=8.6, labelcolor=TEXT,
           bbox_to_anchor=(0.5, -0.04), handlelength=1, columnspacing=1.2)
fig.text(0.125, -0.1, "Differences between healthy people are this large at phylum level, and larger at species level "
         "→ a small case–control study has little power.", fontsize=8.8, color=TEXT2)
header(a, "Metagenomics · the human microbiome", "Every body site is its own ecosystem, and every gut is different",
       dy=40, x1=2.2)
fig.savefig(OUT + "taxonomic_composition.png"); print("ok n1", round(lo, 2), round(hi, 2))
