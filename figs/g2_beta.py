from gallery import *
rng=np.random.default_rng(4)
fig,axes=plt.subplots(1,2,figsize=(11.5,4.6),gridspec_kw={"wspace":0.24})

# (a) β values across CpGs: two modes, few in between
n=60000
beta=np.concatenate([rng.beta(1.4,16,int(n*0.42)),rng.beta(16,1.4,int(n*0.45)),
                     rng.beta(4.5,4.5,int(n*0.13))])
ax=axes[0]
ax.axvspan(0.32,0.68,color=FAINT,lw=0,zorder=0)
ax.hist(beta,bins=110,color=tint(0.45),lw=0,zorder=1)
ax.set_xlabel("β value  (proportion methylated, 0–1)",color=TEXT2)
ax.set_ylabel("number of CpG probes",color=TEXT2)
ax.set_yticks([]); ax.set_xlim(-0.02,1.02)
top=ax.get_ylim()[1]; ax.set_ylim(0,top*1.32)
ax.text(0.03,top*1.28,"unmethylated\nmostly CpG islands, promoters",ha="left",va="top",fontsize=9.5,
        color=TEXT,linespacing=1.35)
ax.text(0.97,top*1.28,"methylated\nmost of the genome",ha="right",va="top",fontsize=9.5,color=TEXT,linespacing=1.35)
ax.text(0.5,top*0.40,"intermediate\nimprinting, X inactivation,\nor a mix of cell types",ha="center",va="center",
        fontsize=9.5,color=TEXT2,linespacing=1.35)
ax.set_title("β values are strongly bimodal",loc="left",fontsize=12.5,fontweight="semibold",color=TEXT,pad=8)

# (b) M = log2(β / (1 − β))
bx=axes[1]
b=np.linspace(0.001,0.999,500)
M=np.log2(b/(1-b))
bx.axhline(0,color=LINE,lw=1); bx.axvline(0.5,color=LINE,lw=1)
bx.plot(b,M,color=CRIMSON,lw=2.2)
bx.set_xlabel("β value",color=TEXT2); bx.set_ylabel("M value  =  log₂( β / (1 − β) )",color=TEXT2)
bx.set_ylim(-7,7); bx.set_xlim(-0.02,1.02)
bx.text(0.04,0.95,"report β,  test on M",transform=bx.transAxes,ha="left",va="top",
        fontsize=11,fontweight="bold",color=CRIMSON)
bx.text(0.96,0.06,"β is bounded: its variance is squeezed at both ends\n→ breaks linear-model assumptions\nM is unbounded, variance roughly constant",
        transform=bx.transAxes,ha="right",va="bottom",fontsize=9.5,color=TEXT2,linespacing=1.5,
        bbox=dict(fc="white",ec="none",pad=2))
bx.set_title("Why analysts convert β to M",loc="left",fontsize=12.5,fontweight="semibold",color=TEXT,pad=8)

for a in axes:
    for s in ("left","bottom"): a.spines[s].set_color(LINE)
    a.tick_params(colors=TEXT2)
header(axes[0],"DNA methylation · β and M values","Report the proportion, test on its logit",dy=24,x1=2.2)
fig.savefig(OUT+"methylation_beta.png"); print("ok g2")
