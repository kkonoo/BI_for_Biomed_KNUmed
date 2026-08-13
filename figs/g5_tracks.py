from style import *
import matplotlib.patches as mp
rng=np.random.default_rng(31)
L=2000
x=np.arange(L)
def bump(c,w,h,sharp=2.0):
    return h*np.exp(-0.5*(np.abs(x-c)/w)**sharp)

prom=[(430,18),(1180,16)]
enh=[(760,42),(1010,38),(1560,50)]

sig={}
sig["H3K4me3"]=sum(bump(c,w,1.0) for c,w in prom)*rng.uniform(0.9,1.0,L)+rng.gamma(1,0.02,L)
sig["H3K27ac"]=(sum(bump(c,w,0.62) for c,w in prom)+sum(bump(c,w,1.0) for c,w in enh))*rng.uniform(0.9,1,L)+rng.gamma(1,0.02,L)
sig["H3K4me1"]=(sum(bump(c,w,0.25) for c,w in prom)+sum(bump(c,w*1.7,0.85) for c,w in enh))*rng.uniform(0.9,1,L)+rng.gamma(1,0.03,L)
sig["H3K27me3"]=bump(1820,190,0.75,1.2)+bump(120,150,0.55,1.2)+rng.gamma(1,0.04,L)
sig["ATAC"]=(sum(bump(c,w*0.55,1.0) for c,w in prom)+sum(bump(c,w*0.5,0.9) for c,w in enh))*rng.uniform(0.9,1,L)+rng.gamma(1,0.015,L)
sig["Input"]=rng.gamma(2.2,0.045,L)+0.06

order=[("ATAC","ATAC-seq","#2E8B57","open chromatin"),
       ("H3K4me3","H3K4me3","#4A77A8","active promoter"),
       ("H3K27ac","H3K27ac","#C9A227","active promoter + enhancer"),
       ("H3K4me1","H3K4me1","#E08A3C","enhancer (primed or active)"),
       ("H3K27me3","H3K27me3","#7B4F9E","polycomb-repressed"),
       ("Input","Input / IgG control","#AAAAAA","background — the control you must have")]

fig,axes=plt.subplots(len(order)+1,1,figsize=(11.5,7.4),sharex=True,
                      gridspec_kw={"height_ratios":[1]*len(order)+[0.85],"hspace":0.0})
for ax,(k,lab,c,note) in zip(axes,order):
    ax.fill_between(x,0,sig[k],color=c,alpha=0.85,lw=0)
    ax.set_ylim(0,1.15); ax.set_xlim(0,L); ax.set_yticks([])
    for sp in ["left","top","right","bottom"]: ax.spines[sp].set_visible(False)
    ax.text(-0.015,0.5,lab,transform=ax.transAxes,ha="right",va="center",
            fontsize=10.5,fontweight="bold",color=c)
    ax.text(1.012,0.5,note,transform=ax.transAxes,ha="left",va="center",fontsize=8.8,color="#777")

gx=axes[-1]
gx.plot([300,1300],[0.66,0.66],color="#5A6570",lw=1.5)
for s,w in [(300,90),(560,60),(840,70),(1130,170)]:
    gx.add_patch(mp.Rectangle((s,0.55),w,0.22,fc=NAVY,ec="none"))
gx.annotate("",xy=(360,0.66),xytext=(300,0.66),arrowprops=dict(arrowstyle="-|>",color=NAVY,lw=1.4))
gx.text(800,0.86,"GENE-A",ha="center",fontsize=9.5,style="italic",color=NAVY)
for c,w in prom: gx.axvspan(c-w*2.2,c+w*2.2,color="#4A77A8",alpha=0.10,lw=0)
for c,w in enh: gx.axvspan(c-w*2.0,c+w*2.0,color="#C9A227",alpha=0.12,lw=0)
gx.text(430,0.30,"promoter",ha="center",fontsize=8.6,color="#4A77A8")
gx.text(1180,0.30,"promoter",ha="center",fontsize=8.6,color="#4A77A8")
for c in [760,1010,1560]:
    gx.text(c,0.30,"enhancer",ha="center",fontsize=8.6,color="#B08D1F")
gx.set_ylim(0,1); gx.set_xlim(0,L); gx.set_yticks([])
for sp in ["left","top","right"]: gx.spines[sp].set_visible(False)
gx.set_xlabel("position (schematic, ~100 kb)")
gx.set_xticks([])

fig.suptitle("Reading epigenomic tracks: the combination of marks defines the element",
             x=0.105,ha="left",fontsize=13,fontweight="bold",y=0.955)
fig.text(0.5,0.02,"H3K4me3 + H3K27ac = active promoter · H3K4me1 + H3K27ac = active enhancer · H3K4me1 alone = primed enhancer · H3K27me3 = repressed",
         ha="center",fontsize=9,color="#666")
fig.savefig(OUT+"browser_tracks.png"); print("ok g5")
