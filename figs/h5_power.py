from style import *
from scipy.stats import norm
fig,axes=plt.subplots(1,2,figsize=(11.5,4.6),gridspec_kw={"wspace":0.28})
n=np.arange(2,31)
ax=axes[0]
for lfc,c,lab in [(0.5,"#C9C9C9","1.4× (log₂FC = 0.5)"),(1.0,BLUE,"2× (log₂FC = 1)"),
                  (2.0,TEAL,"4× (log₂FC = 2)")]:
    cv=0.4
    se=np.sqrt(2)*cv/np.sqrt(n)*1.0
    zcrit=norm.isf(0.05/20000/2)
    pw=norm.sf(zcrit-(lfc*np.log(2)*0.7)/se)
    ax.plot(n,pw,color=c,lw=2.4,label=lab)
ax.axhline(0.8,color="#999",ls="--",lw=1.1)
ax.text(29,0.815,"80% power",ha="right",fontsize=9,color="#777")
ax.axvline(3,color=RED,ls=":",lw=1.6)
ax.text(3.3,0.06,"n = 3",color=RED,fontsize=10,fontweight="bold")
ax.set_xlabel("biological replicates per group"); ax.set_ylabel("power to detect a gene")
ax.set_ylim(0,1.02); ax.set_xlim(2,30)
ax.set_title("Power depends on n far more than on depth",loc="left",fontsize=12)
ax.legend(frameon=False,fontsize=9.2,loc="lower right",title="effect size",title_fontsize=9)

ax=axes[1]
depth=np.array([2,5,10,20,40,80])
reps=np.array([2,3,5,10])
base=np.array([0.30,0.52,0.68,0.76,0.79,0.80])
for r,c,a in zip(reps,[ "#D8D8D8",RED,BLUE,TEAL],[1,1,1,1]):
    scale={2:0.35,3:0.55,5:0.80,10:1.0}[r]
    ax.plot(depth,base*scale,marker="o",ms=5,color=c,lw=2.2,label=f"n = {r}")
ax.set_xscale("log"); ax.set_xticks(depth); ax.set_xticklabels([str(d) for d in depth])
ax.set_xlabel("sequencing depth (million reads per sample)")
ax.set_ylabel("relative sensitivity to detect DE genes")
ax.set_ylim(0,0.9)
ax.set_title("Depth saturates. Replicates do not.",loc="left",fontsize=12)
ax.legend(frameon=False,fontsize=9.2,loc="lower right")
ax.axvspan(20,80,color=OLIVE,alpha=0.13,lw=0)
ax.text(40,0.10,"spending here buys\nalmost nothing",ha="center",fontsize=9,color="#7A6E3C",linespacing=1.4)
ax.legend_.set_bbox_to_anchor((1.0,0.30))
fig.text(0.5,-0.03,"For most bulk RNA-seq questions, 20–30 M reads per sample is enough. Money beyond that is better spent on more samples.",
         ha="center",fontsize=9.3,color="#666")
fig.savefig(OUT+"power_replicates.png"); print("ok h5")
