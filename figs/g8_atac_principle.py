from gallery import *

# ATAC-seq: Tn5 (pre-loaded with adapters) cuts and tags only nucleosome-free DNA
W,Hh=11.6,4.35
fig,ax=plt.subplots(figsize=(11.5,11.5*Hh/W+0.9))
ax.set_xlim(0,W); ax.set_ylim(-1.15,Hh-1.25); ax.set_aspect("equal"); ax.axis("off")
NUC="#D5D9DF"; CW=2.9; DY=1.55

def dna(x0,x1,y,lw=1.6,c=TEXT2): ax.plot([x0,x1],[y,y],color=c,lw=lw,solid_capstyle="butt",zorder=1)
def nuc(x,y,r=0.13): ax.add_patch(mp.Circle((x,y),r,fc=NUC,ec="#AEB4BD",lw=0.8,zorder=2))
def adapter(x,y,d=1):                       # short adapter stub on a DNA end
    ax.add_patch(mp.Rectangle((x if d>0 else x-0.16,y-0.045),0.16,0.09,fc=tint(0.55),ec="none",zorder=4))
def tn5(x,y,s=1.0):                         # Tn5 dimer with two adapter stubs
    for dx in (-0.075*s,0.075*s):
        ax.add_patch(mp.Circle((x+dx,y),0.085*s,fc=CRIMSON,ec="none",zorder=5))
    ax.plot([x-0.12*s,x-0.24*s],[y+0.05*s,y+0.14*s],color=tint(0.55),lw=2.6,zorder=4,solid_capstyle="butt")
    ax.plot([x+0.12*s,x+0.24*s],[y+0.05*s,y+0.14*s],color=tint(0.55),lw=2.6,zorder=4,solid_capstyle="butt")
def chromatin(x0,y,cut=False):
    """closed (packed nucleosomes) on the left, open nucleosome-free stretch, one nucleosome right"""
    dna(x0+0.05,x0+2.65,y)
    for k in range(6): nuc(x0+0.2+k*0.24,y)
    nuc(x0+2.48,y)
    if cut:
        ax.plot([x0+1.72,x0+1.72],[y-0.08,y+0.08],color="white",lw=4,zorder=3)
        ax.plot([x0+2.12,x0+2.12],[y-0.08,y+0.08],color="white",lw=4,zorder=3)
        for xx,d in [(x0+1.72,-1),(x0+1.72,1),(x0+2.12,-1),(x0+2.12,1)]: adapter(xx+0.02*d,y,d)

def head(i,title,where):
    x0=i*CW
    ax.add_patch(mp.Circle((x0+0.2,2.85),0.15,fc=CRIMSON,ec="none"))
    ax.text(x0+0.2,2.85,str(i+1),ha="center",va="center",fontsize=10,fontweight="bold",color="white")
    ax.text(x0+0.45,2.85,title,ha="left",va="center",fontsize=11.5,fontweight="semibold",color=TEXT)
    ax.text(x0+0.06,2.5,where,ha="left",va="center",fontsize=8.5,color=TEXT2,style="italic")
    if i: ax.text(x0-0.12,1.3,"›",ha="center",va="center",fontsize=26,color=LINE)
def caption(i,*lines):
    x0=i*CW
    for k,l in enumerate(lines):
        ax.text(x0+0.06,0.2-k*0.26,l,ha="left",va="top",fontsize=9.2,color=TEXT2)

# 1 intact nucleus
head(0,"Intact nuclei","nucleus kept intact")
ax.add_patch(mp.FancyBboxPatch((0.05,0.55),2.6,1.75,boxstyle="round,pad=0,rounding_size=0.35",
             fc="white",ec=LINE,lw=1,ls=(0,(3,2)),zorder=0))
chromatin(0,1.05)
ax.text(0.8,0.78,"closed",ha="center",va="center",fontsize=9,color=TEXT2)
ax.text(1.92,0.78,"open",ha="center",va="center",fontsize=9,color=CRIMSON,fontweight="semibold")
for xx,yy in [(0.7,1.85),(1.55,1.95),(2.2,1.75)]: tn5(xx,yy)
ax.text(2.42,2.08,"Tn5 + adapters",ha="right",va="bottom",fontsize=9,color=CRIMSON)
caption(0,"Tn5 transposase comes","pre-loaded with","sequencing adapters")

# 2 tagmentation
head(1,"Tagmentation","inside the nucleus")
x0=CW
ax.add_patch(mp.FancyBboxPatch((x0+0.05,0.55),2.6,1.75,boxstyle="round,pad=0,rounding_size=0.35",
             fc="white",ec=LINE,lw=1,ls=(0,(3,2)),zorder=0))
chromatin(x0,1.05,cut=True)
tn5(x0+1.72,1.3,0.8); tn5(x0+2.12,1.3,0.8)
tn5(x0+0.75,1.95,0.8)
ax.annotate("",(x0+0.75,1.3),(x0+0.75,1.8),arrowprops=dict(arrowstyle="-|>",color=LINE,lw=1.1))
ax.text(x0+0.75,1.38,"×",ha="center",va="center",fontsize=14,color=TEXT2,fontweight="bold")
ax.text(x0+0.98,1.6,"no access",ha="left",va="center",fontsize=8.5,color=TEXT2)
ax.text(x0+1.92,0.78,"cut + tag",ha="center",va="center",fontsize=9,color=CRIMSON,fontweight="semibold")
caption(1,"cuts open DNA and pastes","adapters in one step;","packed DNA is not reached")

# 3 fragments -> PCR -> paired-end sequencing
head(2,"Amplify & sequence","in the tube")
x0=2*CW
for k,(a,b) in enumerate([(0.3,1.25),(0.45,1.25),(0.25,1.05)]):
    y=1.95-k*0.28
    dna(x0+a,x0+a+b,y,lw=2)
    adapter(x0+a,y,-1); adapter(x0+a+b,y,1)
ax.text(x0+1.95,1.67,"PCR",ha="left",va="center",fontsize=9,color=TEXT2)
ax.annotate("",(x0+1.35,1.0),(x0+1.35,1.3),arrowprops=dict(arrowstyle="-|>",color=TEXT2,lw=1.1))
dna(x0+0.35,x0+2.35,0.78,lw=2)
adapter(x0+0.35,0.78,-1); adapter(x0+2.35,0.78,1)
ax.annotate("",(x0+0.95,0.62),(x0+0.4,0.62),arrowprops=dict(arrowstyle="-|>",color=CRIMSON,lw=1.6))
ax.annotate("",(x0+1.75,0.62),(x0+2.3,0.62),arrowprops=dict(arrowstyle="-|>",color=CRIMSON,lw=1.6))
ax.text(x0+0.4,0.48,"read 1",ha="left",va="top",fontsize=8.5,color=CRIMSON)
ax.text(x0+2.3,0.48,"read 2",ha="right",va="top",fontsize=8.5,color=CRIMSON)
caption(2,"fragments carry adapters","on both ends → PCR →","paired-end sequencing")

# 4 map reads -> peak over the open region
head(3,"Map & count","on the computer")
x0=3*CW
base=0.62
ax.add_patch(mp.Rectangle((x0+0.08,base-0.02),1.25,1.75,fc=FAINT,ec="none",zorder=0))
ax.text(x0+0.7,2.25,"closed",ha="center",va="center",fontsize=9,color=TEXT2)
ax.text(x0+1.95,2.25,"open",ha="center",va="center",fontsize=9,color=CRIMSON,fontweight="semibold")
rng=np.random.default_rng(3); L=0.2; rows=[]
cent=list(np.clip(rng.normal(x0+1.95,0.3,20),x0+1.4,x0+2.55))+[x0+0.45,x0+1.0,x0+2.55]
for c in sorted(cent):
    a=c-L/2
    for r,end in enumerate(rows):                 # greedy stacking: reads pile up
        if all(a>e+0.03 or a+L<s-0.03 for s,e in end): end.append((a,a+L)); break
    else: rows.append([(a,a+L)]); r=len(rows)-1
    ax.plot([a,a+L],[base+0.08+r*0.11]*2,color=CRIMSON if abs(c-(x0+1.95))<0.5 else tint(0.45),
            lw=2.6,solid_capstyle="butt",zorder=2)
dna(x0+0.08,x0+2.7,base,lw=1.2,c=TEXT2)
ax.text(x0+2.35,base+0.75,"peak",ha="left",va="center",fontsize=9.5,color=CRIMSON,fontweight="semibold")
caption(3,"reads pile up into a peak","over open chromatin;","closed regions stay flat")

ax.plot([0,4*CW-0.2],[-0.68,-0.68],color="#E3E6EB",lw=0.9)
ax.text(0.06,-0.8,"No antibody and no target chosen in advance: one experiment maps promoters, enhancers and other open sites at once. "
        "Schematic.",ha="left",va="top",fontsize=9.5,color=TEXT)
header(ax,"ATAC-seq · assay for transposase-accessible chromatin","Tn5 can only cut where the DNA is open",dy=4)
fig.savefig(OUT+"atac_principle.png"); print("ok g8 atac principle")
