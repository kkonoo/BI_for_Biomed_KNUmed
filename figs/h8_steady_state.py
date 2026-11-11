from gallery import *

# RNA-seq measures a steady-state level: synthesis (tap) and decay (drain) are not separable
W, H = 9.8, 5.0
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, (H - 0.75) / H]); ax.set_xlim(0, W); ax.set_ylim(0, H - 0.75); ax.axis("off")
DARK = "#5E6670"; WATER = tint(0.35)
cases = [("Baseline", 1, 1, "level 100", False),
         ("Gene A: transcription ×2", 2, 1, "level 200", True),
         ("Gene B: half-life ×2", 1, 0.5, "level 200", True)]
TW, TH, BASE = 1.55, 1.75, 1.25                    # tub width, height, bottom y
for k, (name, syn, dec, lvl, hl) in enumerate(cases):
    cx = 1.55 + k * 3.25
    x0, x1 = cx - TW / 2, cx + TW / 2
    level = TH * 0.4 * syn / dec
    ax.add_patch(mp.Rectangle((x0, BASE), TW, level, fc=WATER if hl else "#DDE1E6", ec="none", zorder=1))
    ax.plot([x0, x0, x1, x1], [BASE + TH, BASE, BASE, BASE + TH], color=TEXT2, lw=1.8, zorder=2)
    # tap (transcription)
    ax.plot([x0 - 0.35, x0 + 0.35, x0 + 0.35], [BASE + TH + 0.45, BASE + TH + 0.45, BASE + TH + 0.3],
            color=DARK, lw=3, solid_capstyle="butt", zorder=3)
    ax.add_patch(mp.Rectangle((x0 + 0.35 - 0.035 * syn, BASE + level), 0.07 * syn, TH + 0.3 - level,
                              fc=CRIMSON if syn > 1 else WATER if hl else "#C9CDD3", ec="none", zorder=1))
    # drain (degradation)
    dw = 0.22 * dec
    ax.add_patch(mp.Rectangle((x1 - 0.45 - dw / 2, BASE - 0.02), dw, 0.06, fc="white", ec="none", zorder=3))
    ax.annotate("", (x1 - 0.45, BASE - 0.42), (x1 - 0.45, BASE - 0.02),
                arrowprops=dict(arrowstyle="-|>", color=TEXT2, lw=1.2 + 1.6 * syn,
                                mutation_scale=10))
    # level marker
    ax.plot([x1 + 0.08, x1 + 0.25], [BASE + level, BASE + level], color=TEXT, lw=1)
    ax.text(x1 + 0.3, BASE + level, lvl, ha="left", va="center", fontsize=9.5, fontweight="bold",
            color=CRIMSON if hl else TEXT)
    ax.text(cx, BASE + TH + 0.85, name, ha="center", va="center", fontsize=10.5, fontweight="semibold",
            color=TEXT)
    ax.text(x0 - 0.35, BASE + TH + 0.62, "transcription" + (" (2× faster)" if syn > 1 else ""), ha="left",
            va="bottom", fontsize=8, color=CRIMSON if syn > 1 else TEXT2, fontweight="semibold" if syn > 1 else "normal")
    ax.text(x1 - 0.3, BASE - 0.25, "degradation" + (" (2× slower)" if dec < 1 else ""), ha="left", va="center",
            fontsize=8, color=CRIMSON if dec < 1 else TEXT2, fontweight="semibold" if dec < 1 else "normal")
ax.text(1.55 + TW / 2 + 0.3, BASE + TH * 0.4 - 0.32, "← RNA-seq\n    reads this", ha="left", va="center",
        fontsize=8.3, color=TEXT2, linespacing=1.2)
# verdict
ax.text(W / 2 + 1.6, 0.42, "A and B end at the same level → RNA-seq cannot tell them apart", ha="center",
        va="center", fontsize=9.6, color=CRIMSON, fontweight="semibold")
ax.text(0.1, 0.06, "level = synthesis ÷ decay.  To measure the two separately: nascent RNA (GRO-seq, TT-seq) "
        "or metabolic labelling (SLAM-seq)", ha="left", va="center", fontsize=8.8, color=TEXT2)
header(ax, "Bulk RNA-seq · what it measures", "A count is a water level, not a flow", dy=0, x1=1)
fig.savefig(OUT + "steady_state.png"); print("ok h8 steady state")
