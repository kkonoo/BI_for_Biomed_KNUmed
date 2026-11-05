from gallery import *

# ChIP-seq: crosslink → fragment → immunoprecipitate the mark → sequence; input lane = same chromatin, no antibody
W,Hh=11.6,6.05
fig,ax=plt.subplots(figsize=(11.5,11.5*Hh/W+0.9))
ax.set_xlim(0,W); ax.set_ylim(-1.05,Hh-1.05); ax.set_aspect("equal"); ax.axis("off")
NUC="#D5D9DF"; CW=W/5

def dna(x0,x1,y,lw=1.6,c=TEXT2,z=1): ax.plot([x0,x1],[y,y],color=c,lw=lw,solid_capstyle="butt",zorder=z)
def nuc(x,y,mark=False,r=0.12,alpha=1):
    ax.add_patch(mp.Circle((x,y),r,fc=NUC,ec="#AEB4BD",lw=0.8,zorder=2,alpha=alpha))
    if mark:
        ax.plot([x,x],[y+r,y+r+0.08],color=CRIMSON,lw=1.2,zorder=3,alpha=alpha)
        ax.add_patch(mp.Circle((x,y+r+0.12),0.05,fc=CRIMSON,ec="none",zorder=3,alpha=alpha))
def frag(x,y,mark,alpha=1,w=0.5):
    dna(x-w/2,x+w/2,y,lw=1.5); nuc(x,y,mark,alpha=alpha)
def antibody(x,y,s=1.0,c=CRIMSON):                 # Y with arms pointing down to the target at (x,y)
    ax.plot([x,x],[y+0.12*s,y+0.3*s],color=c,lw=2,solid_capstyle="round",zorder=4)
    ax.plot([x,x-0.09*s],[y+0.12*s,y],color=c,lw=2,solid_capstyle="round",zorder=4)
    ax.plot([x,x+0.09*s],[y+0.12*s,y],color=c,lw=2,solid_capstyle="round",zorder=4)
def head(col,y,n,title,sub):
    x0=col*CW
    ax.add_patch(mp.Circle((x0+0.2,y),0.15,fc=CRIMSON if n.isdigit() else "white",ec=CRIMSON,lw=1.2,zorder=3))
    ax.text(x0+0.2,y,n,ha="center",va="center",fontsize=10,fontweight="bold",
            color="white" if n.isdigit() else CRIMSON,zorder=4)
    ax.text(x0+0.43,y,title,ha="left",va="center",fontsize=11,fontweight="semibold",color=TEXT)
    for k,l in enumerate(sub.split("|")):
        ax.text(x0+0.06,y-0.36-k*0.25,l,ha="left",va="center",fontsize=8.8,color=TEXT2)
def arrow(col,y): ax.text(col*CW-0.1,y,"›",ha="center",va="center",fontsize=24,color=LINE)

YT,YB=4.6,1.45             # lane title heights
ax.text(0.06,YB,"Input control",ha="left",va="center",fontsize=11,color=TEXT,fontweight="semibold")
for k,l in enumerate(["part of the same fragmented","chromatin is set aside","before step 3"]):
    ax.text(0.06,YB-0.36-k*0.25,l,ha="left",va="center",fontsize=8.8,color=TEXT2)

# ① crosslink
head(0,YT,"1","Crosslink","formaldehyde fixes|proteins to DNA")
y=3.2; dna(0.08,2.12,y)
for k,(xx,m) in enumerate([(0.3,False),(0.58,True),(0.86,False),(1.32,True),(1.6,True),(1.95,False)]): nuc(xx,y,m)
ax.add_patch(mp.FancyBboxPatch((1.06,y-0.05),0.14,0.1,boxstyle="round,pad=0.02",fc=TEXT2,ec="none",zorder=2))
ax.text(1.13,y-0.3,"TF",ha="center",va="center",fontsize=8.5,color=TEXT2)
ax.text(1.46,y+0.52,"H3K27ac",ha="center",va="center",fontsize=8.5,color=CRIMSON)
# ② fragment
arrow(1,3.2); head(1,YT,"2","Fragment","sonication or MNase|→ short fragments")
x0=CW
for (dx,dy,m) in [(0.4,3.4,False),(1.1,3.4,True),(1.8,3.4,False),(0.4,2.85,True),(1.1,2.85,False),(1.8,2.85,True)]:
    frag(x0+dx,dy,m,w=0.46)
# ③ immunoprecipitate
arrow(2,3.2); head(2,YT,"3","Immunoprecipitate","an antibody to the mark|pulls down its fragments")
x0=2*CW
ax.add_patch(mp.Circle((x0+1.05,3.62),0.15,fc=TEXT2,ec="none",zorder=4))
ax.text(x0+1.27,3.64,"bead",ha="left",va="center",fontsize=8.5,color=TEXT2)
for dx in (0.6,1.05,1.5):
    frag(x0+dx,2.95,True,w=0.42)
    antibody(x0+dx,2.95+0.29,1.0)
    ax.plot([x0+dx,x0+1.05],[2.95+0.59,3.56],color=CRIMSON,lw=1,zorder=3)
for dx in (0.55,1.45): frag(x0+dx,2.5,False,alpha=0.35,w=0.42)
ax.text(x0+1.0,2.22,"unmarked fragments washed away",ha="center",va="center",fontsize=8.3,color=TEXT2)
# ④ purify & sequence
arrow(3,3.2); head(3,YT,"4","Purify & sequence","reverse crosslinks,|keep the DNA, sequence")
x0=3*CW
for k in range(3):
    yy=3.5-k*0.3; dna(x0+0.35,x0+1.85,yy,lw=2)
ax.annotate("",(x0+0.95,2.62),(x0+0.4,2.62),arrowprops=dict(arrowstyle="-|>",color=CRIMSON,lw=1.5))
ax.annotate("",(x0+1.25,2.62),(x0+1.8,2.62),arrowprops=dict(arrowstyle="-|>",color=CRIMSON,lw=1.5))
ax.text(x0+1.1,2.38,"reads",ha="center",va="center",fontsize=8.5,color=CRIMSON)

# ⑤ map: ChIP and input tracks aligned
arrow(4,3.2); head(4,YT,"5","Map & compare","reads pile up where|the mark sits")
x0=4*CW; xs=np.linspace(x0+0.1,x0+2.2,400); rng=np.random.default_rng(5)
def track(yb,sig,c,lab):
    ax.fill_between(xs,yb,yb+sig,color=c,lw=0,zorder=2); dna(x0+0.1,x0+2.2,yb,lw=0.8,c=LINE)
    ax.text(x0+2.2,yb+0.2,lab,ha="right",va="center",fontsize=8.8,color=c,fontweight="semibold")
bump=lambda m,s,h: h*np.exp(-0.5*((xs-m)/s)**2)
chip=0.03+bump(x0+0.75,0.07,0.6)+bump(x0+1.55,0.09,0.5)+0.02*rng.random(400)
inp=0.05+bump(x0+0.75,0.12,0.08)+bump(x0+1.55,0.15,0.06)+0.03*rng.random(400)
track(2.95,chip,CRIMSON,"ChIP")
track(2.4,inp,"#9AA1AB","input")
ax.text(x0+0.1,2.15,"peak = ChIP enriched over input",ha="left",va="center",fontsize=8.5,color=TEXT,fontweight="semibold")

# input lane: branches off after ②
ax.plot([CW+1.1,CW+1.1,2*CW+0.02],[2.6,YB,YB],color=LINE,lw=1.1,ls=(0,(3,2)),zorder=0)
head(2,YB,"3′","No antibody","same chromatin,|all fragments kept")
x0=2*CW
for (dx,dy,m) in [(0.45,0.3,True),(1.05,0.3,False),(1.65,0.3,True),(0.75,-0.2,False),(1.35,-0.2,False)]:
    frag(x0+dx,dy,m,w=0.42)
arrow(3,0.15); head(3,YB,"4′","Purify & sequence","same as ④")
x0=3*CW
for k in range(3): dna(x0+0.35,x0+1.85,0.4-k*0.3,lw=2,c="#9AA1AB")
ax.text(4*CW+0.1,YB,"Background only:",ha="left",va="center",fontsize=9.5,color=TEXT,fontweight="semibold")
ax.text(4*CW+0.1,YB-0.33,"open regions shear more easily,",ha="left",va="center",fontsize=8.8,color=TEXT2)
ax.text(4*CW+0.1,YB-0.6,"copy-number gains add reads",ha="left",va="center",fontsize=8.8,color=TEXT2)
ax.text(4*CW+0.1,YB-0.95,"IgG control: same steps with a",ha="left",va="center",fontsize=8.8,color=TEXT2)
ax.text(4*CW+0.1,YB-1.22,"non-specific antibody",ha="left",va="center",fontsize=8.8,color=TEXT2)

ax.plot([0,W],[-0.55,-0.55],color="#E3E6EB",lw=0.9)
ax.text(0.02,-0.7,"Without a matched input or IgG control, enrichment cannot be told apart from fragmentation and copy-number artefacts. Schematic.",
        ha="left",va="top",fontsize=9.5,color=TEXT)
header(ax,"ChIP-seq · chromatin immunoprecipitation","An antibody pulls out the DNA that carries one mark",dy=4)
fig.savefig(OUT+"chip_principle.png"); print("ok g9 chip principle")
