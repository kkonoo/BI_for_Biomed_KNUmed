from gallery import *
from math import lgamma, exp

# ORA background: the same overlap tested against expressed genes vs the whole annotation
N_ALL, N_EXPR, K_PATH, K_PATH_EXPR, N_DE, OVER = 20000, 12000, 200, 190, 600, 14

def lchoose(n, k): return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)
def p_greater(k, N, K, n):
    """One-sided Fisher's exact (hypergeometric upper tail) P(X >= k)."""
    return sum(exp(lchoose(K, i) + lchoose(N - K, n - i) - lchoose(N, n)) for i in range(k, min(K, n) + 1))

p_ok = p_greater(OVER, N_EXPR, K_PATH_EXPR, N_DE)
p_bad = p_greater(OVER, N_ALL, K_PATH, N_DE)
e_ok = N_DE * K_PATH_EXPR / N_EXPR
e_bad = N_DE * K_PATH / N_ALL

fig = plt.figure(figsize=(9.8, 5.0))
gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.05], wspace=0.12, left=0.01, right=0.99, top=0.8, bottom=0.12)
DARK = "#5E6670"; GREYF = "#B8BEC7"

def panel(ax, t, sub):
    ax.set_title(t, loc="left", fontsize=12.5, fontweight="semibold", color=TEXT, pad=22)
    ax.text(0, 1.02, sub, transform=ax.transAxes, ha="left", va="bottom", fontsize=9.3, color=TEXT2)

# (a) the gene universe ---------------------------------------------------------------------------------
a = fig.add_subplot(gs[0, 0]); a.set_xlim(0, 10); a.set_ylim(0, 7.2); a.set_aspect("equal"); a.axis("off")
panel(a, "(a) Which genes could be on the list?", "the DE list can only come from expressed genes")
a.add_patch(mp.FancyBboxPatch((0.1, 0.1), 9.8, 7.0, boxstyle="round,pad=0,rounding_size=0.25",
                              fc=FAINT, ec=LINE, lw=1))
a.text(0.4, 6.75, "all annotated genes", ha="left", va="center", fontsize=10, fontweight="semibold", color=TEXT)
a.text(0.4, 6.35, f"{N_ALL:,}", ha="left", va="center", fontsize=9.5, color=TEXT2)
a.add_patch(mp.FancyBboxPatch((0.4, 0.4), 6.6, 5.6, boxstyle="round,pad=0,rounding_size=0.25",
                              fc=tint(0.1), ec=tint(0.55), lw=1.3))
a.text(0.7, 5.65, "expressed in this tissue", ha="left", va="center", fontsize=10, fontweight="semibold",
       color=CRIMSON)
a.text(0.7, 5.25, f"{N_EXPR:,}  ← correct background", ha="left", va="center", fontsize=9.5, color=CRIMSON)
a.text(8.45, 3.3, f"{N_ALL - N_EXPR:,}\nnot expressed:\nnever\ndetectable", ha="center", va="center",
       fontsize=9, color=TEXT2, linespacing=1.35)
# DE list and pathway circles
a.add_patch(mp.Circle((2.9, 2.55), 1.55, fc="white", ec=DARK, lw=1.3, alpha=0.95))
a.text(2.35, 2.6, f"DE list\n{N_DE}", ha="center", va="center", fontsize=9.5, color=TEXT, fontweight="semibold",
       linespacing=1.3)
a.add_patch(mp.Circle((4.75, 2.55), 1.05, fc=tint(0.35), ec=CRIMSON, lw=1.3, alpha=0.75))
a.text(5.15, 2.55, f"pathway\n{K_PATH}", ha="center", va="center", fontsize=9.5, color=TEXT, fontweight="semibold",
       linespacing=1.3)
a.text(4.05, 2.55, f"{OVER}", ha="center", va="center", fontsize=11, color=CRIMSON, fontweight="bold")
a.annotate("overlap", (4.05, 2.2), (4.05, 0.85), ha="center", va="center", fontsize=9, color=CRIMSON,
           arrowprops=dict(arrowstyle="-", color=CRIMSON, lw=0.9))
a.text(5.35, 3.95, f"{K_PATH_EXPR} of {K_PATH} expressed", ha="center", va="center", fontsize=8.6, color=TEXT2)

# (b) the same overlap, two backgrounds ------------------------------------------------------------------
b = fig.add_subplot(gs[0, 1]); b.set_xlim(0, 10); b.set_ylim(0, 7.2); b.axis("off")
panel(b, "(b) Same overlap, two backgrounds", f"observed overlap = {OVER} genes in both cases")
def case(y, ok, lab, bg, exp_, p):
    col = CRIMSON if ok else DARK
    b.text(0.0, y + 1.6, ("✓  " if ok else "✕  ") + lab, ha="left", va="center", fontsize=10.5,
           fontweight="semibold", color=col)
    b.text(0.0, y + 1.15, f"background = {bg:,} genes", ha="left", va="center", fontsize=9, color=TEXT2)
    sc = 0.34
    for k, (name, v, fc) in enumerate((("expected", exp_, GREYF), ("observed", OVER, tint(0.55) if ok else DARK))):
        yy = y + 0.45 - k * 0.55
        b.text(0.0, yy, name, ha="left", va="center", fontsize=9, color=TEXT2)
        b.add_patch(mp.Rectangle((1.6, yy - 0.17), v * sc, 0.34, fc=fc, ec="none"))
        b.text(1.7 + v * sc, yy, f"{v:.1f}" if name == "expected" else f"{v}", ha="left", va="center",
               fontsize=9, color=TEXT)
    b.text(10, y + 0.2, f"{OVER / exp_:.1f}× enriched\np = {p:.3f}", ha="right", va="center", fontsize=10.5,
           fontweight="bold", color=col, linespacing=1.35)
case(4.0, True, "expressed genes", N_EXPR, e_ok, p_ok)
case(0.55, False, "whole annotation", N_ALL, e_bad, p_bad)
b.text(10, 3.55, "not significant", ha="right", va="center", fontsize=9, color=CRIMSON)
b.text(10, 0.1, "looks significant — but isn't", ha="right", va="center", fontsize=9, color=DARK)

fig.text(0.01, 0.06, f"Adding the {N_ALL - N_EXPR:,} genes that could never be detected lowers the expected overlap "
         "→ every enrichment p-value is inflated.", ha="left", va="center", fontsize=9.3, color=TEXT2)
a.apply_aspect(); pa, pb = a.get_position(), b.get_position()
b.set_position([pb.x0, pa.y0, pb.width, pa.height])
header(a, "Enrichment analysis · ORA background", "Background = the genes that could have been detected",
       dy=40, x1=2.08)
fig.savefig(OUT + "ora_background.png"); print("ok h10 background", p_ok, p_bad, e_ok, e_bad)
