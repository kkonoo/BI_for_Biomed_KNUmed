from style import *
rng=np.random.default_rng(4)
fig,axes=plt.subplots(1,2,figsize=(11.5,4.6),gridspec_kw={"wspace":0.26})
n=60000
beta=np.concatenate([rng.beta(1.4,16,int(n*0.42)),rng.beta(16,1.4,int(n*0.45)),
                     rng.beta(4.5,4.5,int(n*0.13))])
ax=axes[0]
ax.hist(beta,bins=110,color=TEAL,alpha=0.75,lw=0)
ax.set_xlabel("β value  (proportion methylated, 0–1)"); ax.set_ylabel("number of CpG probes")
ax.set_title("β values are strongly bimodal",loc="left",fontsize=12)
ax.set_yticks([])
for x,t,c in [(0.06,"unmethylated\n(mostly CpG islands\nand promoters)","#2E8B57"),
              (0.94,"methylated\n(most of the genome)","#7B4F9E")]:
    ax.text(x,ax.get_ylim()[1]*0.80,t,ha="center",fontsize=9.3,color=c,linespacing=1.4)
ax.axvspan(0.32,0.68,color=OLIVE,alpha=0.16,lw=0)
ax.text(0.5,ax.get_ylim()[1]*0.34,"intermediate:\nimprinting, X-inactivation,\nor mixed cell types",
        ha="center",fontsize=9,color="#7A6E3C",linespacing=1.4)

ax=axes[1]
b=np.linspace(0.001,0.999,500)
M=np.log2(b/(1-b))
ax.plot(b,M,color=PURPLE,lw=2.2)
ax.axhline(0,color="#CCC",lw=1); ax.axvline(0.5,color="#CCC",lw=1)
ax.set_xlabel("β value"); ax.set_ylabel("M value  =  log₂( β / (1−β) )")
ax.set_title("Why analysts convert β to M",loc="left",fontsize=12)
ax.set_ylim(-7,7)
ax.text(0.5,-5.4,"β is bounded and its variance is\ncompressed at both ends — bad for\nlinear models.  M is unbounded and\nroughly homoscedastic.",
        ha="center",fontsize=9.3,color="#444",linespacing=1.5,
        bbox=dict(boxstyle="round,pad=0.45",fc="#FAFAFA",ec="#DDD"))
ax.text(0.04,0.94,"report β,  test on M",transform=ax.transAxes,ha="left",va="top",
        fontsize=11,fontweight="bold",color=PURPLE)
fig.savefig(OUT+"methylation_beta.png"); print("ok g2")
