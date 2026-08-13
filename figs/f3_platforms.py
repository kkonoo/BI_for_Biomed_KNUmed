from style import *
import matplotlib.patches as mp
rng=np.random.default_rng(7)
fig,axes=plt.subplots(4,1,figsize=(10,6.0),sharex=True,
                      gridspec_kw={"height_ratios":[0.55,1,1,1],"hspace":0.42})
L=1000
exons=[(80,55),(230,40),(330,70),(520,45),(690,60),(860,50)]

ga=axes[0]
ga.plot([0,L],[0.5,0.5],color="#9A9A9A",lw=1.6,zorder=1)
for s,w in exons:
    ga.add_patch(mp.Rectangle((s,0.22),w,0.56,fc=NAVY,ec="none",zorder=2))
ga.set_ylim(0,1); ga.set_xlim(0,L); ga.axis("off")
ga.text(-8,0.5,"gene",ha="right",va="center",fontsize=10.5,fontweight="bold")
ga.text(L+8,0.5,"exons (dark) · introns (line)",ha="left",va="center",fontsize=8.5,color="#777")

x=np.arange(L)
inex=np.zeros(L,bool)
for s,w in exons: inex[s:s+w]=True

def track(ax,y,color,title,note):
    ax.fill_between(x,0,y,color=color,alpha=0.55,lw=0)
    ax.plot(x,y,color=color,lw=0.8)
    ax.set_ylim(0,105); ax.set_xlim(0,L)
    ax.set_yticks([0,50,100]); ax.set_ylabel("depth",fontsize=9)
    ax.set_title(title,loc="left",fontsize=11.5,color=color)
    ax.text(0.995,0.86,note,transform=ax.transAxes,ha="right",va="top",
            fontsize=9,color="#666",style="italic",
            bbox=dict(boxstyle="round,pad=0.28",fc="white",ec="none",alpha=0.88))

wgs=30+rng.normal(0,3.2,L); wgs=np.clip(wgs,0,None)
track(axes[1],wgs,BLUE,"WGS — whole genome sequencing",
      "even coverage everywhere · finds SVs and non-coding variants · most expensive")

wes=np.where(inex,95+rng.normal(0,16,L),rng.normal(1.2,0.8,L)); wes=np.clip(wes,0,None)
for i in range(1,L):
    if not inex[i] and inex[i-1]: wes[i:i+18]=np.linspace(wes[i-1],1.5,min(18,L-i))
track(axes[2],wes,ORANGE,"WES — whole exome sequencing",
      "deep in exons, blind outside · ~2% of the genome · uneven capture efficiency")

arr=np.zeros(L)
sites=np.sort(rng.choice(L,42,replace=False))
track(axes[3],arr,TEAL,"Genotyping array",
      "only pre-selected sites are measured · no new variants · cheapest by far")
axes[3].vlines(sites,0,70,color=TEAL,lw=1.8)
axes[3].scatter(sites,[70]*len(sites),s=14,color=TEAL,zorder=3)
axes[3].set_ylabel("genotyped",fontsize=9); axes[3].set_yticks([])
axes[3].set_xlabel("position along the gene (bp, schematic)")

for a in axes[1:]:
    for s,w in exons:
        a.axvspan(s,s+w,color=NAVY,alpha=0.05,lw=0)
fig.suptitle("What each platform actually measures",x=0.09,ha="left",fontsize=13,fontweight="bold",y=0.985)
fig.savefig(OUT+"platform_comparison.png"); print("ok f3")
