from gallery import *
import os
from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch
from matplotlib.font_manager import FontProperties
from matplotlib.transforms import Affine2D

rng=np.random.default_rng(7)
fig,(ax,bx)=plt.subplots(1,2,figsize=(11.5,5.4),gridspec_kw={"width_ratios":[1.3,1],"wspace":0.12})

# (a) peaks vs matched background: each line = one region, crimson ticks = motif hits
N=22
def stack(x0,hits,title,sub):
    rows=rng.permutation(N)[:hits]
    for i in range(N):
        y=6.5-i*0.3; L=rng.uniform(3.0,3.6); s=x0+rng.uniform(0,3.6-L)
        ax.plot([s,s+L],[y,y],color=LINE,lw=2.2,solid_capstyle="round")
        if i in rows:
            for _ in range(rng.choice([1,1,2])):
                xm=rng.uniform(s+0.2,s+L-0.2)
                ax.plot([xm,xm],[y-0.11,y+0.11],color=CRIMSON,lw=2.6)
    ax.text(x0,7.9,title,ha="left",va="bottom",fontsize=11.5,fontweight="semibold",color=TEXT)
    ax.text(x0,7.82,sub,ha="left",va="top",fontsize=9.5,color=TEXT2,linespacing=1.35)

stack(0.0,8,"Peaks","1,000 ATAC / ChIP peaks")
stack(4.6,2,"Background","1,000 regions matched for\nGC content and openness")
for x,pc,frac,c in [(0.0,"38%","380 of 1,000 contain the motif",CRIMSON),
                    (4.6,"9%","90 of 1,000 contain the motif",TEXT2)]:
    ax.text(x,-0.1,pc,ha="left",va="top",fontsize=24,fontweight="bold",color=c)
    ax.text(x,-1.1,frac,ha="left",va="top",fontsize=9.5,color=TEXT2)
ax.plot([0,8.2],[-1.75,-1.75],color="#E3E6EB",lw=0.9)
ax.text(0,-2.0,"Enrichment test (hypergeometric / Fisher's exact)",ha="left",va="top",
        fontsize=10.5,fontweight="semibold",color=TEXT)
ax.text(0,-2.5,"4.2× more often than in the background  →  p ≈ 10⁻⁵⁵",ha="left",va="top",
        fontsize=10.5,color=CRIMSON,fontweight="semibold")
ax.add_patch(mp.Rectangle((2.55,7.97),0.07,0.3,fc=CRIMSON,ec="none"))
ax.text(2.72,8.12,"= motif hit",ha="left",va="center",fontsize=9.5,color=TEXT2)
ax.set_xlim(-0.1,8.3); ax.set_ylim(-3.15,8.6); ax.axis("off")
header(ax,"Downstream analysis · TF motif enrichment",
       "Which TF motifs are over-represented in my peaks?",dy=26,x1=1.875)
ax.set_title("(a) Count motif hits in peaks and in a background set",loc="left",
             fontsize=12.5,fontweight="semibold",color=TEXT,pad=8)

# (b) sequence logo of the enriched motif (AP-1, TGA[C/G]TCA)
FP=FontProperties(fname=os.path.expanduser("~/.fonts/Pretendard-Bold.otf"))
COL={"A":CRIMSON,"T":tint(0.55),"C":"#8A8F98","G":"#3F4248"}
PWM=[dict(T=.88,C=.04,A=.04,G=.04), dict(G=.90,A=.06,C=.02,T=.02), dict(A=.92,C=.03,G=.03,T=.02),
     dict(C=.48,G=.46,A=.03,T=.03), dict(T=.90,C=.05,A=.03,G=.02), dict(C=.88,A=.06,G=.04,T=.02),
     dict(A=.86,G=.06,C=.04,T=.04)]
def letter(ch,x,y,w,h):
    tp=TextPath((0,0),ch,size=1,prop=FP); b=tp.get_extents()
    tr=Affine2D().translate(-b.x0,-b.y0).scale(w/b.width,h/b.height).translate(x,y)
    bx.add_patch(PathPatch(tr.transform_path(tp),fc=COL[ch],ec="none"))
for i,p in enumerate(PWM):
    ic=2+sum(v*np.log2(v) for v in p.values())
    y=0
    for ch,v in sorted(p.items(),key=lambda kv:kv[1]):
        h=v*ic
        if h>0.02: letter(ch,i+0.06,y,0.88,h)
        y+=h
bx.plot([0,7],[0,0],color=LINE,lw=1)
bx.plot([-0.15,-0.15],[0,2],color=LINE,lw=1)
for v in (0,1,2):
    bx.plot([-0.22,-0.15],[v,v],color=LINE,lw=1)
    bx.text(-0.3,v,str(v),ha="right",va="center",fontsize=9.5,color=TEXT2)
bx.text(-0.75,1,"bits",ha="center",va="center",rotation=90,fontsize=9.5,color=TEXT2)
for i in range(7): bx.text(i+0.5,-0.1,str(i+1),ha="center",va="top",fontsize=9.5,color=TEXT2)
bx.set_xlim(-0.9,7.1); bx.set_ylim(-1.75,2.35); bx.axis("off")
bx.set_title("(b) The enriched motif: AP-1  (TGA[C/G]TCA)",loc="left",fontsize=12.5,
             fontweight="semibold",color=TEXT,pad=8)
bx.text(3.5,1.25,"C or G",ha="center",va="bottom",fontsize=9,color=TEXT2)
for k,(t1,t2) in enumerate([("AP-1 = FOS, JUN, ATF … all share this motif →",
                              "the motif alone cannot say which one is bound"),
                             ("The background set decides the answer:",
                              "random genome regions inflate every enrichment")]):
    y=-0.55-k*0.6
    bx.text(-0.15,y,t1,ha="left",va="top",fontsize=9.5,color=TEXT if k==0 else TEXT2)
    bx.text(-0.15,y-0.21,t2,ha="left",va="top",fontsize=9.5,color=TEXT2)
fig.savefig(OUT+"motif_enrichment.png"); print("ok g21 motif")
