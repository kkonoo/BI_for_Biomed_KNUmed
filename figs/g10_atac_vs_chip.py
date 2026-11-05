from gallery import *

# Why ChIP does not need accessibility: ATAC works in intact nuclei, ChIP fragments chromatin first
fig=plt.figure(figsize=(11.5,6.9))
ax=fig.add_axes([0,0.33,1,0.6]); bx=fig.add_axes([0,0.0,1,0.29])
W=11.6; ax.set_xlim(0,W); ax.set_ylim(-0.2,3.6); ax.set_aspect("equal"); ax.axis("off")
NUC="#D5D9DF"; HET="#B9BFC8"; DARK="#5E6670"

def dna(x0,x1,y,lw=1.5,c=TEXT2,z=1): ax.plot([x0,x1],[y,y],color=c,lw=lw,solid_capstyle="butt",zorder=z)
def nuc(x,y,r=0.12,fc=NUC,flag=None):
    ax.add_patch(mp.Circle((x,y),r,fc=fc,ec="#A6ADB7",lw=0.8,zorder=2))
    if flag:
        ax.plot([x,x],[y+r,y+r+0.07],color=flag,lw=1.1,zorder=3)
        ax.add_patch(mp.Circle((x,y+r+0.11),0.045,fc=flag,ec="none",zorder=3))
def hetero(cx,cy):                          # compact heterochromatin: tightly packed nucleosomes
    pts=[(cx+dx*0.2,cy+dy*0.19) for dy in (-1,0,1) for dx in (-1.5,-0.5,0.5,1.5) if not (abs(dy)==1 and abs(dx)==1.5)]
    ax.add_patch(mp.FancyBboxPatch((cx-0.48,cy-0.33),0.96,0.66,boxstyle="round,pad=0.02,rounding_size=0.3",
                 fc="#EEF0F3",ec="none",zorder=1))
    for x,y in pts: nuc(x,y,r=0.1,fc=HET)
def tn5(x,y,s=0.8):
    for dx in (-0.075*s,0.075*s): ax.add_patch(mp.Circle((x+dx,y),0.085*s,fc=CRIMSON,ec="none",zorder=5))
    ax.plot([x-0.12*s,x-0.24*s],[y+0.05*s,y+0.14*s],color=tint(0.55),lw=2.4,zorder=4)
    ax.plot([x+0.12*s,x+0.24*s],[y+0.05*s,y+0.14*s],color=tint(0.55),lw=2.4,zorder=4)
def antibody(x,y,c=CRIMSON):
    ax.plot([x,x],[y+0.12,y+0.28],color=c,lw=2,solid_capstyle="round",zorder=4)
    ax.plot([x,x-0.09],[y+0.12,y],color=c,lw=2,solid_capstyle="round",zorder=4)
    ax.plot([x,x+0.09],[y+0.12,y],color=c,lw=2,solid_capstyle="round",zorder=4)
def nucleus(x0,y0,w,h):
    ax.add_patch(mp.FancyBboxPatch((x0,y0),w,h,boxstyle="round,pad=0,rounding_size=0.45",
                 fc="white",ec=LINE,lw=1.1,ls=(0,(3,2)),zorder=0))

def panel_title(x,t,sub):
    ax.text(x,3.4,t,ha="left",va="center",fontsize=12.5,fontweight="semibold",color=TEXT)
    ax.text(x,3.08,sub,ha="left",va="center",fontsize=9.5,color=TEXT2)

# ATAC (left)
panel_title(0.1,"ATAC-seq","Tn5 must reach the DNA inside an intact nucleus")
nucleus(0.15,0.25,4.9,2.45)
hetero(1.25,1.25); ax.text(1.25,0.62,"heterochromatin",ha="center",va="center",fontsize=9,color=TEXT2)
dna(2.35,4.85,1.25)
for xx in (2.55,2.8,4.4,4.65): nuc(xx,1.25)
ax.text(3.6,0.95,"open element",ha="center",va="center",fontsize=9,color=CRIMSON,fontweight="semibold")
tn5(3.45,1.45); tn5(3.8,1.45)
tn5(1.25,2.3); ax.annotate("",(1.25,1.7),(1.25,2.15),arrowprops=dict(arrowstyle="-|>",color=LINE,lw=1.1))
ax.text(1.25,1.76,"×",ha="center",va="center",fontsize=15,fontweight="bold",color=TEXT2,zorder=6)
ax.text(1.5,2.0,"cannot enter",ha="left",va="center",fontsize=8.8,color=TEXT2)
ax.text(3.62,2.0,"cuts here",ha="center",va="center",fontsize=8.8,color=CRIMSON)
ax.text(0.15,0.0,"→ closed chromatin gives no ATAC signal",ha="left",va="center",fontsize=10,color=TEXT,fontweight="semibold")

ax.plot([5.75,5.75],[0.0,3.5],color="#E3E6EB",lw=0.9)

# ChIP (right)
panel_title(6.1,"ChIP-seq","chromatin is crosslinked and fragmented first")
nucleus(6.15,0.55,1.6,1.85)
hetero(6.95,1.75)
dna(6.3,7.6,1.0); nuc(6.6,1.0); nuc(7.3,1.0)
ax.annotate("",(8.4,1.48),(7.85,1.48),arrowprops=dict(arrowstyle="-|>",color=TEXT2,lw=1.2))
ax.text(8.12,1.72,"fragment",ha="center",va="center",fontsize=8.8,color=TEXT2)
# tube of fragments in solution
ax.add_patch(mp.FancyBboxPatch((8.5,0.35),2.95,2.25,boxstyle="round,pad=0,rounding_size=0.3",
             fc="#FAFAFB",ec=LINE,lw=1,zorder=0))
ax.text(10.97,2.42,"in solution",ha="center",va="center",fontsize=8.5,color=TEXT2,style="italic")
for x,y,f in [(8.95,1.9,True),(9.75,2.05,True),(10.55,1.75,False),(9.1,1.0,False),(9.95,1.15,True),(10.8,0.85,True)]:
    dna(x-0.22,x+0.22,y,lw=1.4); nuc(x,y,r=0.1,fc=HET if f else NUC,flag=DARK if f else None)
    if f: antibody(x,y+0.27)
ax.text(8.6,0.52,"antibody binds H3K9me3 fragments",ha="left",va="center",fontsize=8.6,color=TEXT2)
ax.text(6.1,0.0,"→ marks on closed chromatin can be mapped",ha="left",va="center",fontsize=10,color=TEXT,fontweight="semibold")

header(ax,"Accessibility","ATAC needs open DNA; ChIP does not",dy=6)

# bottom: same locus, aligned tracks
xs=np.linspace(0,100,800); rng=np.random.default_rng(8)
bump=lambda m,s,h: h*np.exp(-0.5*((xs-m)/s)**2)
atac=0.02+bump(62,1.6,1.0)+bump(84,1.2,0.45)+0.02*rng.random(800)
k9=0.06+0.75/(1+np.exp(-(xs-14)/1.5))*(1/(1+np.exp((xs-44)/1.5)))*(0.85+0.3*rng.random(800))
X0,X1=14,96
def tr(yb,sig,c,lab,sub):
    x=X0+(X1-X0)*xs/100
    bx.fill_between(x,yb,yb+0.75*sig,color=c,lw=0,zorder=2)
    bx.plot([X0,X1],[yb,yb],color=LINE,lw=0.8)
    bx.text(X0-1,yb+0.3,lab,ha="right",va="center",fontsize=10.5,fontweight="semibold",color=c)
    bx.text(X0-1,yb+0.05,sub,ha="right",va="center",fontsize=8.8,color=TEXT2)
bx.add_patch(mp.Rectangle((X0+0.14*(X1-X0),-0.15),0.3*(X1-X0),2.15,fc=FAINT,ec="none",zorder=0))
tr(1.05,atac,CRIMSON,"ATAC","accessibility")
tr(0.1,k9,DARK,"ChIP H3K9me3","heterochromatin mark")
bx.text(X0+0.29*(X1-X0),2.12,"heterochromatin",ha="center",va="center",fontsize=9,color=TEXT2)
bx.text(X0+0.62*(X1-X0),2.12,"open element",ha="center",va="center",fontsize=9,color=CRIMSON)
bx.set_xlim(0,100); bx.set_ylim(-0.25,2.3); bx.axis("off")
bx.text(X0,-0.2,"Same locus, schematic. The two assays complement each other: ATAC finds open elements, "
        "ChIP maps a mark wherever it sits.",ha="left",va="top",fontsize=9.5,color=TEXT)
fig.savefig(OUT+"atac_vs_chip.png"); print("ok g10 atac vs chip")
