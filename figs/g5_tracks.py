from gallery import *
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

# ATAC = the dark emphasis; active marks = crimson tints; repressive = dark grey; input = light grey
order=[("ATAC","ATAC-seq",CRIMSON,"open chromatin"),
       ("H3K4me3","H3K4me3",tint(0.75),"active promoter"),
       ("H3K27ac","H3K27ac",tint(0.55),"active promoter or enhancer"),
       ("H3K4me1","H3K4me1",tint(0.38),"enhancer (primed or active)"),
       ("H3K27me3","H3K27me3","#5F6673","Polycomb-repressed"),
       ("Input","Input / IgG","#C9CED6","background (control)")]

fig,axes=plt.subplots(len(order)+1,1,figsize=(11.5,7.2),sharex=True,
                      gridspec_kw={"height_ratios":[1]*len(order)+[0.95],"hspace":0.08})
for ax,(k,lab,c,note) in zip(axes,order):
    ax.fill_between(x,0,sig[k],color=c,lw=0)
    ax.axhline(0,color=LINE,lw=0.8)
    ax.set_ylim(0,1.15); ax.set_xlim(0,L); ax.set_yticks([])
    for sp in ["left","top","right","bottom"]: ax.spines[sp].set_visible(False)
    ax.text(-0.015,0.42,lab,transform=ax.transAxes,ha="right",va="center",
            fontsize=10.5,fontweight="semibold",color=TEXT)
    ax.text(1.012,0.42,note,transform=ax.transAxes,ha="left",va="center",fontsize=9.5,color=TEXT2)

# gene model and element annotation
gx=axes[-1]
for c,w in prom: gx.axvspan(c-w*2.2,c+w*2.2,ymin=0.0,ymax=1,color=FAINT,lw=0)
for c,w in enh: gx.axvspan(c-w*2.0,c+w*2.0,ymin=0.0,ymax=1,color=FAINT,lw=0)
def gene(tss,exons,end,name):
    gx.plot([tss,end],[0.66,0.66],color=TEXT2,lw=1.3)
    for s0,w in exons: gx.add_patch(mp.Rectangle((s0,0.55),w,0.22,fc=TEXT,ec="none"))
    gx.plot([tss,tss],[0.77,0.95],color=TEXT,lw=1.2)
    gx.annotate("",xy=(tss+70,0.95),xytext=(tss,0.95),arrowprops=dict(arrowstyle="-|>",color=TEXT,lw=1.2))
    gx.text(tss+85,0.95,name,ha="left",va="center",fontsize=9.5,style="italic",color=TEXT)
gene(422,[(422,55),(600,45),(890,50),(1060,70)],1130,"GENE-A")
gene(1172,[(1172,45),(1300,40),(1430,45)],1475,"GENE-B")
for c,_ in prom: gx.text(c,0.22,"promoter",ha="center",fontsize=9,color=TEXT,fontweight="semibold")
for c,_ in enh:  gx.text(c,0.22,"enhancer",ha="center",fontsize=9,color=TEXT2)
gx.set_ylim(0,1.1); gx.set_xlim(0,L); gx.set_yticks([]); gx.set_xticks([])
for sp in ["left","top","right"]: gx.spines[sp].set_visible(False)
gx.spines["bottom"].set_color(LINE)
gx.set_xlabel("position  (schematic, ~100 kb)",color=TEXT2)

# shade element columns lightly across all tracks
for ax in axes[:-1]:
    for c,w in prom: ax.axvspan(c-w*2.2,c+w*2.2,color=FAINT,lw=0,zorder=0)
    for c,w in enh: ax.axvspan(c-w*2.0,c+w*2.0,color=FAINT,lw=0,zorder=0)

fig.text(0.125,0.035,"H3K4me3 + H3K27ac = active promoter  ·  H3K4me1 + H3K27ac = active enhancer  ·  "
         "H3K4me1 alone = primed enhancer  ·  H3K27me3 = repressed",ha="left",fontsize=9.5,color=TEXT)
header(axes[0],"Reading epigenomic tracks","The combination of marks, not one track, defines the element",dy=8)
fig.savefig(OUT+"browser_tracks.png"); print("ok g5")
