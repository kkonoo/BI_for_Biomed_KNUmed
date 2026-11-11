from gallery import *

# Protein inference: peptides map to proteins many-to-many; what is reported depends on which peptides were seen
W, H = 9.8, 4.6
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, (H - 0.8) / H]); ax.set_xlim(0, W); ax.set_ylim(0, H - 0.8); ax.axis("off")
DARK = "#5E6670"; GREYF = "#B8BEC7"
prot = {"A": 0.9, "B": 2.1, "C": 3.3}
pep = {f"p{i}": 0.35 + (i - 1) * 0.66 for i in range(1, 7)}
links = {"A": ["p1", "p2", "p3"], "B": ["p3", "p4", "p5"], "C": ["p5", "p6"]}
owners = {p: [k for k, v in links.items() if p in v] for p in pep}
YP, YQ = 2.8, 1.0

def pill(x, y, label, fc, ec, tc, w=0.5, h=0.3, fs=8.5):
    ax.add_patch(mp.FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0,rounding_size=0.12",
                                   fc=fc, ec=ec, lw=1))
    ax.text(x, y, label, ha="center", va="center", fontsize=fs, color=tc, fontweight="semibold")

# (a) the map
for pr, ps in links.items():
    for p in ps:
        shared = len(owners[p]) > 1
        ax.plot([prot[pr], pep[p]], [YP - 0.2, YQ + 0.15], color=CRIMSON if shared else GREYF,
                lw=1.4 if shared else 1.2, ls=(0, (3, 2)) if shared else "-", zorder=1)
for pr, x in prot.items():
    pill(x, YP, f"protein {pr}", DARK, DARK, "white", w=0.95, h=0.4, fs=9.5)
for p, x in pep.items():
    shared = len(owners[p]) > 1
    pill(x, YQ, p, "white", CRIMSON if shared else LINE, CRIMSON if shared else TEXT)
ax.text(0.1, 0.45, "grey = unique to one protein", ha="left", va="center", fontsize=8.8, color=TEXT2)
ax.text(0.1, 0.18, "red dashed = shared (paralogues, isoforms)", ha="left", va="center", fontsize=8.8, color=CRIMSON)
ax.text(0.1, H - 1.25, "(a) Peptides map to proteins many-to-many", ha="left", va="center", fontsize=12.3,
        fontweight="semibold", color=TEXT)

# (b) three cases
X0 = 4.55
ax.text(X0, H - 1.25, "(b) What gets reported depends on what you saw", ha="left", va="center", fontsize=12.3,
        fontweight="semibold", color=TEXT)
ax.text(X0, H - 1.55, "filled = detected peptides", ha="left", va="center", fontsize=9, color=TEXT2)
cases = [(["p1", "p2", "p3"], "protein A", "two unique peptides → solid", TEXT),
         (["p3"], "group {A, B}", "only a shared peptide → A or B?", CRIMSON),
         (["p5", "p6"], "protein C", "one unique peptide → weak", DARK)]
for k, (seen, res, why, col) in enumerate(cases):
    y = 2.5 - k * 0.86
    for i, p in enumerate(pep):
        xx = X0 + 0.2 + i * 0.36
        on = p in seen
        shared = len(owners[p]) > 1
        ax.add_patch(mp.FancyBboxPatch((xx - 0.15, y - 0.11), 0.3, 0.22, boxstyle="round,pad=0,rounding_size=0.08",
                                       fc=(CRIMSON if shared else DARK) if on else "white",
                                       ec=(CRIMSON if shared else DARK) if on else LINE, lw=0.9))
        ax.text(xx, y - 0.24, p, ha="center", va="center", fontsize=6.8, color=TEXT2)
    ax.annotate("", (X0 + 2.5, y), (X0 + 2.25, y), arrowprops=dict(arrowstyle="-|>", color=GREYF, lw=1,
                                                                     mutation_scale=8))
    ax.text(X0 + 2.65, y + 0.1, res, ha="left", va="center", fontsize=10.5, fontweight="bold", color=col)
    ax.text(X0 + 2.65, y - 0.17, why, ha="left", va="center", fontsize=8.6, color=TEXT2)
ax.text(X0, 0.12, "Check: how many unique peptides support the protein you care about?", ha="left", va="center",
        fontsize=9, color=CRIMSON, fontweight="semibold")
header(ax, "Proteomics · protein inference", "One shared peptide is not a protein identification", x=0.1 / W,
       dy=0, x1=1 - 0.1 / W)
fig.savefig(OUT + "protein_inference.png"); print("ok m4 inference")
