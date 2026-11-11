from gallery import *
rng = np.random.default_rng(17)

# Reverse causation (simulated). Two stories give the same snapshot at diagnosis; only samples taken before onset
# tell them apart. (a) snapshot; (b) the microbe falls before onset; (c) the microbe falls after onset.
t = np.linspace(-24, 12, 73)
n = 12
sig = lambda x, c: 1 / (1 + np.exp(-(x - c) / 1.8))
def people(centre):
    return np.array([-1.6 * sig(t, centre + rng.normal(0, 1.5)) + rng.normal(0, 0.18, len(t)).cumsum() * 0.12
                     + rng.normal(0, 0.15) for _ in range(n)])
cause, conseq = people(-9), people(1.5)
ctrl = np.array([rng.normal(0, 0.18, len(t)).cumsum() * 0.12 + rng.normal(0, 0.15) for _ in range(n)])
t_s = 6                                                      # sampled at diagnosis, some months after onset
i_s = np.argmin(np.abs(t - t_s))

fig, axes = plt.subplots(1, 3, figsize=(9.8, 4.0), gridspec_kw={"wspace": 0.28, "width_ratios": [0.65, 1, 1]})
DARK = "#5E6670"; GREYF = "#B8BEC7"

def style(ax, tt, sub):
    ax.set_title(tt, loc="left", fontsize=12.3, fontweight="semibold", color=TEXT, pad=22)
    ax.text(0, 1.03, sub, transform=ax.transAxes, fontsize=9, color=TEXT2)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"): ax.spines[sp].set_color(LINE)
    ax.tick_params(colors=TEXT2, labelsize=8.3); ax.set_ylim(-2.6, 1.0)

a = axes[0]
for x, v, col in ((0, ctrl[:, i_s], GREYF), (1, np.r_[cause[:, i_s]], CRIMSON)):
    a.scatter(x + rng.uniform(-0.15, 0.15, n), v, s=16, color=col, lw=0)
    a.plot([x - 0.25, x + 0.25], [np.median(v)] * 2, color=TEXT, lw=2)
a.set_xticks([0, 1]); a.set_xticklabels(["controls", "patients"], fontsize=9, color=TEXT)
a.tick_params(axis="x", length=0); a.set_xlim(-0.6, 1.6)
a.set_ylabel("Faecalibacterium (log, relative to start)", color=TEXT2, fontsize=8.6)
style(a, "(a) One snapshot", "sampled at diagnosis")

for ax, data, tt, sub in ((axes[1], cause, "(b) Falls before onset", "compatible with a cause"),
                          (axes[2], conseq, "(c) Falls after onset", "a consequence of disease")):
    ax.axvspan(0, 12, color=FAINT, lw=0, zorder=0)
    for row in data:
        ax.plot(t, row, color=CRIMSON, lw=0.9, alpha=0.55)
    ax.axvline(0, color=TEXT, lw=0.8)
    ax.text(-0.5, 0.85, "disease\nonset", ha="right", va="top", fontsize=8.3, color=TEXT, linespacing=1.05)
    ax.axvline(t_s, color=TEXT2, lw=0.7, ls=(0, (3, 2)))
    ax.text(t_s + 0.4, -2.45, "snapshot", ha="left", va="bottom", fontsize=8, color=TEXT2)
    ax.set_xlim(-24, 12); ax.set_xticks([-24, -12, 0, 12])
    ax.set_xlabel("months from onset", color=TEXT2, fontsize=9)
    style(ax, tt, sub)
    ax.set_yticklabels([])
header(a, "Metagenomics · reverse causation", "A snapshot cannot tell cause from consequence; a time series can",
       dy=40, x1=4.25)
fig.text(0.125, -0.07, "Simulated. Both stories give the same picture at diagnosis (a). Only samples taken before onset "
         "separate (b) from (c).", fontsize=8.8, color=TEXT2)
fig.savefig(OUT + "longitudinal.png"); print("ok n17")
