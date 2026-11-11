from gallery import *

# Hi-C QC on one page: valid pairs, cis/trans, long-range cis, P(s)
fig = plt.figure(figsize=(9.8, 8.8))
gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.16, left=0.01, right=0.99, top=0.85, bottom=0.05)
DARK = "#5E6670"; GREYF = "#B8BEC7"

def panel(ax, t, sub):
    ax.set_title(t, loc="left", fontsize=12.5, fontweight="semibold", color=TEXT, pad=22)
    ax.text(0, 1.02, sub, transform=ax.transAxes, ha="left", va="bottom", fontsize=9.3, color=TEXT2)

# (a) valid pairs vs artefacts --------------------------------------------------------------------------
a = fig.add_subplot(gs[0, 0]); a.set_xlim(0, 10); a.set_ylim(-0.3, 4.6); a.axis("off")
panel(a, "(a) Valid pairs", "keep only reads that join two different fragments")
def frag(x0, x1, y, c, lw=4): a.plot([x0, x1], [y, y], color=c, lw=lw, solid_capstyle="butt")
def read(x, y, d, c):
    a.annotate("", (x + 0.75 * d, y), (x, y), arrowprops=dict(arrowstyle="-|>", color=c, lw=1.3, mutation_scale=8))
def junction(x, y):
    a.add_patch(mp.Circle((x, y), 0.09, fc="white", ec=CRIMSON, lw=1.1, zorder=4))
# valid
y = 3.75
frag(0.4, 2.3, y, CRIMSON); frag(2.3, 4.2, y, DARK); junction(2.3, y)
read(0.6, y + 0.3, 1, CRIMSON); read(4.0, y + 0.3, -1, DARK)
a.text(4.6, y + 0.12, "✓ valid pair", ha="left", va="center", fontsize=10, fontweight="semibold", color=CRIMSON)
a.text(4.6, y - 0.28, "fragment i + fragment j, ligated → a real contact", ha="left", va="center",
       fontsize=8.8, color=TEXT2)
# dangling end
y = 2.25
frag(0.4, 4.2, y, GREYF)
a.plot([4.2, 4.2], [y - 0.2, y + 0.2], color=DARK, lw=1)
read(0.6, y + 0.3, 1, DARK); read(4.0, y + 0.3, -1, DARK)
a.text(4.36, y + 0.42, "cut end", ha="left", va="center", fontsize=8, color=TEXT2)
a.text(4.6, y + 0.0, "✕ dangling end", ha="left", va="center", fontsize=10, fontweight="semibold", color=DARK)
a.text(4.6, y - 0.4, "one fragment, never ligated", ha="left", va="center", fontsize=8.8, color=TEXT2)
# self-circle
y = 0.75
a.add_patch(mp.Circle((2.3, y), 0.55, fc="none", ec=GREYF, lw=4))
junction(2.3, y + 0.55)
a.annotate("", (1.55, y + 0.25), (1.7, y - 0.35), arrowprops=dict(arrowstyle="-|>", color=DARK, lw=1.3,
                                                                  mutation_scale=8))
a.annotate("", (3.05, y + 0.25), (2.9, y - 0.35), arrowprops=dict(arrowstyle="-|>", color=DARK, lw=1.3,
                                                                  mutation_scale=8))
a.text(4.6, y + 0.2, "✕ self-ligation", ha="left", va="center", fontsize=10, fontweight="semibold", color=DARK)
a.text(4.6, y - 0.2, "a fragment ligated to itself", ha="left", va="center", fontsize=8.8, color=TEXT2)

# (b) cis vs trans --------------------------------------------------------------------------------------
b = fig.add_subplot(gs[0, 1]); b.set_xlim(0, 10); b.set_ylim(-0.3, 4.6); b.axis("off")
panel(b, "(b) Cis / trans ratio", "contacts within one chromosome vs between chromosomes")
for y0, lab in ((3.7, "chr A"), (2.5, "chr B")):
    b.add_patch(mp.FancyBboxPatch((0.6, y0 - 0.12), 4.0, 0.24, boxstyle="round,pad=0,rounding_size=0.12",
                                  fc=FAINT, ec=LINE, lw=0.8))
    b.text(0.45, y0, lab, ha="right", va="center", fontsize=9, color=TEXT2)
for x0, x1 in ((1.0, 2.0), (1.6, 3.9), (2.7, 3.4)):
    xs = np.linspace(x0, x1, 50)
    b.plot(xs, 3.82 + 0.4 * np.sin(np.pi * (xs - x0) / (x1 - x0)), color=CRIMSON, lw=1.4)
b.text(4.8, 4.05, "cis", ha="left", va="center", fontsize=10, color=CRIMSON, fontweight="semibold")
for x0, x1 in ((1.3, 3.2), (3.6, 2.0)):
    b.plot([x0, x1], [3.58, 2.62], color=DARK, lw=1.2, ls=(0, (3, 2)))
b.text(4.8, 3.1, "trans", ha="left", va="center", fontsize=10, color=DARK, fontweight="semibold")
def bar(y, cis, lab, good):
    b.add_patch(mp.Rectangle((0.6, y), 6.0 * cis, 0.36, fc=CRIMSON if good else tint(0.45), ec="none"))
    b.add_patch(mp.Rectangle((0.6 + 6.0 * cis, y), 6.0 * (1 - cis), 0.36, fc=GREYF, ec="none"))
    b.text(0.6, y + 0.5, lab, ha="left", va="bottom", fontsize=9.6, fontweight="semibold",
           color=TEXT if good else TEXT2)
bar(1.15, 0.85, "good: mostly cis", True)
bar(0.0, 0.5, "poor: lots of trans → random ligation", False)
b.text(6.8, 1.33, "cis | trans", ha="left", va="center", fontsize=8.8, color=TEXT2)

# (c) long-range cis -----------------------------------------------------------------------------------
c = fig.add_subplot(gs[1, 0])
panel(c, "(c) Long-range cis fraction", "cis contacts farther than 20 kb carry the 3D information")
lx = np.linspace(2, 8, 600)                      # log10 distance (100 bp – 100 Mb)
dens = 0.55 * np.exp(-0.5 * ((lx - 2.9) / 0.35) ** 2) + 0.45 * np.exp(-0.5 * ((lx - 5.3) / 1.0) ** 2)
cut = np.log10(20000)
c.fill_between(lx, dens, where=lx < cut, color=GREYF, lw=0)
c.fill_between(lx, dens, where=lx >= cut, color=tint(0.45), lw=0)
c.plot(lx, dens, color=TEXT2, lw=1)
c.axvline(cut, color=TEXT, lw=1, ls=(0, (3, 2)))
c.text(cut + 0.08, 0.6, "20 kb", ha="left", va="top", fontsize=9, color=TEXT)
c.text(2.05, 0.79, "short-range: mostly\nligation by-products", ha="left", va="top", fontsize=8.8,
       color=TEXT2, linespacing=1.3)
c.text(6.2, 0.36, "long-range cis\n→ loops, TADs,\ncompartments", ha="left", va="center", fontsize=8.8,
       color=CRIMSON, fontweight="semibold", linespacing=1.3)
c.set_xticks([2, 3, 4, 5, 6, 7, 8]); c.set_xticklabels(["100 bp", "1 kb", "10 kb", "100 kb", "1 Mb", "10 Mb", "100 Mb"])
c.set_xlim(2, 8); c.set_ylim(0, 0.82); c.set_yticks([])
c.set_xlabel("distance between the two ends (cis pairs)", color=TEXT2, fontsize=9)
c.set_ylabel("read pairs", color=TEXT2, fontsize=9)

# (d) P(s) --------------------------------------------------------------------------------------------
d = fig.add_subplot(gs[1, 1])
panel(d, "(d) Distance decay P(s)", "contact probability falls smoothly with distance")
s = np.logspace(4, 8, 200)
good = (s / 1e4) ** -1.0
poor = (s / 1e4) ** -1.0 + 3e-3
d.loglog(s, good, color=CRIMSON, lw=1.8)
d.loglog(s, poor, color=TEXT2, lw=1.5, ls=(0, (4, 2)))
d.text(2e5, 0.12, "good: a straight line\n(slope ≈ −1 in the Mb range)", ha="left", va="center", fontsize=8.8,
       color=CRIMSON, fontweight="semibold", linespacing=1.3)
d.text(2.5e6, 7e-3, "poor: flattens far away\n→ random-ligation floor", ha="left", va="bottom", fontsize=8.8,
       color=TEXT2, linespacing=1.3)
d.set_xlim(1e4, 1e8); d.set_ylim(5e-5, 2)
d.set_xticks([1e4, 1e5, 1e6, 1e7, 1e8]); d.set_xticklabels(["10 kb", "100 kb", "1 Mb", "10 Mb", "100 Mb"])
d.set_yticks([]); d.minorticks_off()
d.set_xlabel("genomic distance s (log)", color=TEXT2, fontsize=9)
d.set_ylabel("contact probability (log)", color=TEXT2, fontsize=9)
for ax in (c, d):
    for sp in ("left", "bottom"): ax.spines[sp].set_color(LINE)
    ax.tick_params(colors=TEXT2, labelsize=8.5)

header(a, "Hi-C QC", "Four checks: did proximity ligation really happen?", dy=40, x1=2.03)
fig.savefig(OUT + "hic_qc.png"); print("ok g26 hic qc")
