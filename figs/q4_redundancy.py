from gallery import *
import os
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "assets", "img", "w7") + "/"

# Functional redundancy (illustrative numbers): two healthy guts share none of their top genera,
# yet their pathway profiles are almost the same. Butyrate producers (crimson) differ by name, not by job.
DARK = "#5E6670"; G1 = "#8E959E"; G2 = "#B8BEC7"; G3 = "#D5D9DF"; G4 = "#E9ECEF"
taxa = {"Person 1": [("Bacteroides", 38, DARK), ("Faecalibacterium", 18, CRIMSON), ("Alistipes", 12, G1),
                     ("Akkermansia", 8, G2), ("Bifidobacterium", 6, G3), ("other genera", 18, G4)],
        "Person 2": [("Prevotella", 42, DARK), ("Roseburia", 10, CRIMSON), ("E. rectale", 9, tint(0.6)),
                     ("Ruminococcus", 11, G1), ("Dialister", 6, G2), ("other genera", 22, G4)]}
paths = {"Person 1": [("carbohydrate breakdown", 31, DARK), ("amino acid synthesis", 22, G1),
                      ("butyrate production", 9, CRIMSON), ("cell wall and membrane", 14, G2),
                      ("vitamin B synthesis", 7, G3), ("other pathways", 17, G4)],
         "Person 2": [("carbohydrate breakdown", 30, DARK), ("amino acid synthesis", 23, G1),
                      ("butyrate production", 9, CRIMSON), ("cell wall and membrane", 13, G2),
                      ("vitamin B synthesis", 7, G3), ("other pathways", 18, G4)]}

fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.9), gridspec_kw={"wspace": 0.55})
def stacks(a, data, side_labels):
    xs = [0, 1.25]
    for x, (person, segs) in zip(xs, data.items()):
        bottom = 0
        for name, v, col in segs:
            a.bar(x, v, bottom=bottom, width=0.55, color=col, edgecolor="white", lw=1.2)
            dark = col in (DARK, CRIMSON, G1) or col == tint(0.6)
            a.text(x, bottom + v / 2, f"{v}%", ha="center", va="center", fontsize=8.4,
                   color="white" if dark else TEXT, fontweight="semibold" if col == CRIMSON else "normal")
            if (x == xs[0] and side_labels == "left") or x == xs[1]:
                lx, ha = (x - 0.36, "right") if x == xs[0] else (x + 0.36, "left")
                a.text(lx, bottom + v / 2, name, ha=ha, va="center", fontsize=8.6,
                       color=CRIMSON if col in (CRIMSON, tint(0.6)) else TEXT,
                       style="italic" if a is axes[0] and "other" not in name else "normal",
                       fontweight="semibold" if col in (CRIMSON, tint(0.6)) else "normal")
            bottom += v
    a.set_xticks(xs); a.set_xticklabels(list(data.keys()), fontsize=9.5, color=TEXT)
    a.set_ylim(0, 100); a.set_xlim(-1.3, 2.6)
    a.set_yticks([]); a.tick_params(axis="x", length=0)
    for sp in ("top", "right", "left"): a.spines[sp].set_visible(False)
    a.spines["bottom"].set_color(LINE)

a = axes[0]; stacks(a, taxa, "left")
a.set_title("(a) Who is there", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
a.text(0, 1.03, "share of reads by genus — no top genus in common", transform=a.transAxes, fontsize=9, color=TEXT2)
b = axes[1]; stacks(b, paths, "none")
b.set_xlim(-0.55, 3.25)
b.set_title("(b) What they can do", loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
b.text(0, 1.03, "share of reads by pathway (HUMAnN-style)", transform=b.transAxes, fontsize=9, color=TEXT2)
fig.patches.append(mp.FancyArrowPatch((0.425, 0.47), (0.535, 0.47), transform=fig.transFigure,
                                      arrowstyle="-|>", mutation_scale=14, color="#9AA1AB", lw=1.3))
fig.text(0.48, 0.50, "read the genes,\nsum into pathways", ha="center", va="bottom", fontsize=8.6, color=TEXT2,
         linespacing=1.15)
fig.text(0.5, -0.02, "Red = butyrate producers. Person 1 relies on Faecalibacterium, person 2 on Roseburia and "
         "E. rectale — different names, the same job. Illustrative numbers.",
         ha="center", fontsize=8.8, color=TEXT2)
header(a, "Metagenomics · functional profiling", "Different taxa, same functions: functional redundancy",
       dy=40, x1=2.55)
fig.savefig(OUT + "functional_redundancy.png"); print("ok q4 redundancy")
