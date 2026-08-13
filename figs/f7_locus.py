from style import *
import matplotlib.patches as mp
from matplotlib.lines import Line2D
rng=np.random.default_rng(21)
fig,(ax,gx)=plt.subplots(2,1,figsize=(10.5,6.0),sharex=True,
                         gridspec_kw={"height_ratios":[3.1,1.0],"hspace":0.12})
lo,hi=44.0,44.9
lead=44.46
n=1400
pos=np.sort(rng.uniform(lo,hi,n))
d=np.abs(pos-lead)
r2=np.exp(-((d/0.055)**1.7))*rng.uniform(0.75,1.0,n)
r2=np.clip(r2,0,1)
logp=r2*13.5+rng.gamma(1.0,0.85,n)
logp=np.clip(logp,0.02,None)
i=np.argmin(np.abs(pos-lead)); logp[i]=logp.max()+0.7; r2[i]=1.0

bins=[(0.8,1.01,"#C1544B","0.8 – 1.0"),(0.6,0.8,"#E08A3C","0.6 – 0.8"),
      (0.4,0.6,"#C9A227","0.4 – 0.6"),(0.2,0.4,"#5F9EA0","0.2 – 0.4"),
      (-0.01,0.2,"#5C79A8","< 0.2")]
for a,b,c,lab in bins[::-1]:
    m=(r2>a)&(r2<=b)&(np.arange(n)!=i)
    ax.scatter(pos[m],logp[m],s=26,color=c,alpha=0.85,lw=0.3,edgecolor="white")
ax.scatter([pos[i]],[logp[i]],s=170,marker="D",color="#7B2D8E",ec="white",lw=1.4,zorder=6)
ax.annotate("lead SNP\nrs1234567",xy=(pos[i],logp[i]),xytext=(pos[i]-0.115,logp[i]+0.6),
            fontsize=10,fontweight="bold",color="#7B2D8E",ha="center")
ax.axhline(-np.log10(5e-8),color=RED,lw=1.2,ls="--")
ax.text(hi-0.005,-np.log10(5e-8)+0.25,"p = 5×10⁻⁸",ha="right",fontsize=9,color=RED)
ax.set_ylabel(r"$-\log_{10}(p)$"); ax.set_ylim(0,17.4); ax.set_xlim(lo,hi)
ax.set_title("Regional association plot — an association is a region, not a point",loc="left")
ax.legend(handles=[Line2D([],[],marker="o",ls="",color=c,label=lab,markersize=7) for _,_,c,lab in bins],
          title="LD  $r^2$  with lead SNP",frameon=False,fontsize=8.8,title_fontsize=9,
          loc="upper left",handletextpad=0.3,borderpad=0.2)
ax.text(0.985,0.93,"every point here is correlated with\nevery other one — you cannot tell\nwhich is causal from this figure alone",
        transform=ax.transAxes,ha="right",va="top",fontsize=9.3,color="#444",linespacing=1.5,
        bbox=dict(boxstyle="round,pad=0.45",fc="#FAFAFA",ec="#DDDDDD"))

genes=[("GENE1",44.06,44.19,1),("LNC-AS1",44.28,44.34,-1),("GENE2",44.41,44.62,1),
       ("GENE3",44.70,44.79,-1),("GENE4",44.80,44.88,1)]
for j,(nm,s,e,strand) in enumerate(genes):
    y=0.72 if j%2==0 else 0.28
    gx.add_patch(mp.FancyArrow(s if strand>0 else e,y,(e-s)*strand,0,width=0.055,
                 head_width=0.13,head_length=0.018,length_includes_head=True,
                 fc=NAVY if nm!="GENE2" else "#7B2D8E",ec="none"))
    gx.text((s+e)/2,y+0.15,nm,ha="center",fontsize=8.6,style="italic",
            color=NAVY if nm!="GENE2" else "#7B2D8E")
gx.set_ylim(0,1); gx.set_xlim(lo,hi); gx.set_yticks([])
gx.set_xlabel("Position on chromosome 6 (Mb)")
for sp in ["left","top","right"]: gx.spines[sp].set_visible(False)
fig.text(0.5,-0.01,"The nearest gene to the lead SNP is not necessarily the causal gene — see §4.3.",
         ha="center",fontsize=9,color="#666")
fig.savefig(OUT+"locus_ld.png"); print("ok f7")
