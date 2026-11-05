from gallery import *
rng=np.random.default_rng(16)

fig,ax=plt.subplots(figsize=(11.5,6.9))
x=np.linspace(0,100,4000)
g=lambda c,s: np.exp(-0.5*((x-c)/s)**2)
def plateau(a,b,edge=1.2): return 1/(1+np.exp(-(x-a)/edge))*(1/(1+np.exp((x-b)/edge)))

ENH,TSS,TES=12,26,55; SIL=(62,79); REP=(84,100)
DARK="#5A5F66"; DFILL="#A3A8B0"
tracks=[  # name, signal, active?, (label x, label y), shape → meaning
 ("H3K4me3", 1.0*g(TSS+0.6,0.8), True, (TSS+2.2,0.78),"sharp peak at the TSS","active promoter"),
 ("H3K27ac", 0.85*g(ENH,1.1)+0.75*g(TSS-1.6,0.9)+0.45*g(TSS+2.4,0.9), True,(TSS+4.6,0.62),
  "peaks at promoter and enhancer","active promoter or enhancer"),
 ("H3K4me1", 0.75*g(ENH,2.4)+0.38*g(TSS-4,1.8)+0.38*g(TSS+5,1.9), True,(ENH+4.6,0.62),
  "broader peaks: enhancer, promoter shoulders","enhancer (primed or active)"),
 ("H3K36me3", np.clip((x-(TSS+2))/(TES-TSS-2),0,1)**1.3*0.8*plateau(TSS+2,TES,0.8), True,(TES-11,0.82),
  "rises along the gene body toward 3′","transcribed gene body"),
 ("H3K27me3", 0.62*plateau(*SIL,1.0), False,(SIL[0]+1,0.8),"broad domain over a silent gene","Polycomb-repressed"),
 ("H3K9me3", 0.7*plateau(*REP,1.5), False,(REP[0]+2.5,0.86),"very broad domain over repeats",
  "constitutive heterochromatin"),
]
RH=1.15                                       # row height
for i,(name,sig,act,(lx,ly),shape,meaning) in enumerate(tracks):
    y0=-i*RH
    s=sig+0.035+0.03*rng.random(x.size)*(1+0.6*(sig>0.05))
    if not act:
        n=np.convolve(rng.standard_normal(x.size),np.ones(30)/30,mode="same")*0.9
        s=s*(1+n*(sig>0.05))
    ax.fill_between(x,y0,y0+0.8*s,color=tint(0.5) if act else DFILL,lw=0,zorder=2)
    ax.plot(x,y0+0.8*s,color=CRIMSON if act else DARK,lw=0.7,zorder=3)
    ax.plot([0,100],[y0,y0],color=LINE,lw=0.8,zorder=1)
    ax.text(lx,y0+ly,name,ha="left",va="center",fontsize=10.5,fontweight="bold",
            color=CRIMSON if act else DARK,zorder=4,bbox=dict(fc="white",ec="none",pad=0.6,alpha=0.85))
    ax.text(104,y0+0.5,shape,ha="left",va="center",fontsize=9.5,color=TEXT2)
    ax.text(104,y0+0.2,"→ "+meaning,ha="left",va="center",fontsize=10,fontweight="semibold",
            color=CRIMSON if act else TEXT)

# element columns + locus model
yb=-len(tracks)*RH+0.35
for (a,b) in [(ENH-3,ENH+3),(TSS-3,TSS+3),SIL,REP]:
    ax.axvspan(a,b,ymin=0,ymax=1,color=FAINT,zorder=0)
ax.plot([0,100],[yb,yb],color=TEXT2,lw=1.1,zorder=2)
ax.add_patch(mp.Rectangle((ENH-1.8,yb-0.12),3.6,0.24,fc=tint(0.35),ec="none",zorder=3))
for a,b in [(TSS,TSS+2),(33,34.4),(40,41.2),(47,48.2),(TES-2.4,TES)]:
    ax.add_patch(mp.Rectangle((a,yb-0.17),b-a,0.34,fc=TEXT,ec="none",zorder=3))
ax.plot([TSS,TSS],[yb+0.17,yb+0.45],color=TEXT,lw=1.1)
ax.annotate("",(TSS+3.6,yb+0.45),(TSS,yb+0.45),arrowprops=dict(arrowstyle="-|>",color=TEXT,lw=1.1,mutation_scale=10))
for a,b in [(SIL[0]+1,SIL[0]+2.6),(68,69.2),(73,74.2),(SIL[1]-2.5,SIL[1]-1)]:
    ax.add_patch(mp.Rectangle((a,yb-0.17),b-a,0.34,fc=DFILL,ec="none",zorder=3))
for k in range(6):
    ax.add_patch(mp.Rectangle((REP[0]+1+k*2.6,yb-0.1),1.8,0.2,fc="white",ec=DARK,lw=0.8,hatch="////",zorder=3))
labs=[(ENH,"enhancer",TEXT),((TSS+TES)/2,"active gene (TSS → gene body)",TEXT),
      (sum(SIL)/2,"silent gene",TEXT),(sum(REP)/2,"repeats",TEXT)]
for xm,t,c in labs:
    ax.text(xm,yb-0.42,t,ha="center",va="top",fontsize=10,fontweight="semibold",color=c)
ax.text(sum(REP)/2,yb-0.78,"(heterochromatin)",ha="center",va="top",fontsize=9,color=TEXT2)
ax.text(104,yb,"schematic locus, ~100 kb",ha="left",va="center",fontsize=9,color=TEXT2)
ax.text(106,1.05,"active marks",ha="left",va="bottom",fontsize=9.5,color=CRIMSON,fontweight="semibold")
ax.text(123,1.05,"repressive marks",ha="left",va="bottom",fontsize=9.5,color=DARK,fontweight="semibold")
ax.add_patch(mp.Rectangle((104,1.1),1.4,0.22,fc=tint(0.5),ec="none"))
ax.add_patch(mp.Rectangle((121,1.1),1.4,0.22,fc=DFILL,ec="none"))
ax.set_xlim(-2,141); ax.set_ylim(yb-1.05,1.45); ax.axis("off")
header(ax,"Histone modifications","Each mark has its own place and its own shape",dy=0)
fig.savefig(OUT+"histone_marks.png"); print("ok g16 histone marks")
