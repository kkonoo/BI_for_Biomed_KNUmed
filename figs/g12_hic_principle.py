from gallery import *
from matplotlib.path import Path

# Hi-C: crosslink → digest → biotin fill-in → proximity ligation → pull-down + paired-end seq → contact matrix
W,Hh=11.6,7.4
fig,ax=plt.subplots(figsize=(11.5,11.5*Hh/W+0.9))
ax.set_xlim(0,W); ax.set_ylim(-0.85,Hh-0.85); ax.set_aspect("equal"); ax.axis("off")
DARK="#5E6670"; PROT="#B9BFC8"; CW=W/3
CI,CJ=CRIMSON,DARK                       # locus i, locus j

def seg(x0,x1,y,c,lw=3.2,z=2): ax.plot([x0,x1],[y,y],color=c,lw=lw,solid_capstyle="butt",zorder=z)
def biotin(x,y):
    ax.add_patch(mp.Circle((x,y),0.075,fc="white",ec=CRIMSON,lw=1.2,zorder=5))
    ax.text(x,y,"B",ha="center",va="center",fontsize=6.5,color=CRIMSON,fontweight="bold",zorder=6)
def prot(x,y,w=0.32,h=0.5): ax.add_patch(mp.FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle="round,pad=0,rounding_size=0.12",
                                                          fc=PROT,ec="none",zorder=1))
def head(r,c,n,title,sub):
    x0=c*CW; y=6.25 if r==0 else 2.85
    ax.add_patch(mp.Circle((x0+0.2,y),0.15,fc=CRIMSON,ec="none"))
    ax.text(x0+0.2,y,str(n),ha="center",va="center",fontsize=10,fontweight="bold",color="white")
    ax.text(x0+0.43,y,title,ha="left",va="center",fontsize=11.5,fontweight="semibold",color=TEXT)
    for k,l in enumerate(sub.split("|")):
        ax.text(x0+0.06,y-0.36-k*0.25,l,ha="left",va="center",fontsize=9,color=TEXT2)

# ① crosslink: a looping fibre brings i and j together
head(0,0,1,"Crosslink","formaldehyde fixes loci that|touch in 3D")
# hairpin loop: i and j run side by side, the loop between them sticks out to the right
P=Path([(2.4,4.85),(3.0,4.85),(3.55,5.55),(3.45,4.7),(3.4,3.85),(3.0,4.5),(2.4,4.5)],
       [Path.MOVETO]+[Path.CURVE4]*6)
ax.add_patch(mp.PathPatch(P,fc="none",ec=LINE,lw=2.2,zorder=1))
ax.plot([0.3,1.5],[5.2,4.85],color=LINE,lw=2.2,zorder=1); ax.plot([0.3,1.5],[4.15,4.5],color=LINE,lw=2.2,zorder=1)
prot(1.95,4.675,0.32,0.55)
seg(1.5,2.4,4.85,CI); seg(1.5,2.4,4.5,CJ)
ax.text(1.38,4.95,"i",ha="right",va="center",fontsize=11,color=CI,fontweight="bold")
ax.text(1.38,4.38,"j",ha="right",va="center",fontsize=11,color=CJ,fontweight="bold")
# linear view
ax.plot([0.3,3.6],[3.65,3.65],color=LINE,lw=1.6)
seg(0.75,1.05,3.65,CI,lw=5); seg(3.0,3.3,3.65,CJ,lw=5)
ax.text(0.9,3.85,"i",ha="center",va="center",fontsize=9.5,color=CI,fontweight="bold")
ax.text(3.15,3.85,"j",ha="center",va="center",fontsize=9.5,color=CJ,fontweight="bold")
ax.text(1.95,3.42,"far apart on the linear genome",ha="center",va="center",fontsize=8.5,color=TEXT2)

# ② digest
ax.text(CW-0.1,4.8,"›",ha="center",va="center",fontsize=24,color=LINE)
head(0,1,2,"Digest","a restriction enzyme cuts|(DpnII / MboI: GATC)")
x0=CW; prot(x0+1.9,4.75,0.4,0.62)
seg(x0+0.5,x0+1.25,5.0,CI); seg(x0+1.38,x0+1.85,5.0,CI); seg(x0+1.98,x0+2.45,5.0,CI); seg(x0+2.58,x0+3.2,5.0,CI)
seg(x0+0.5,x0+1.05,4.5,CJ); seg(x0+1.18,x0+1.85,4.5,CJ); seg(x0+1.98,x0+2.6,4.5,CJ); seg(x0+2.73,x0+3.2,4.5,CJ)
for xx in (1.31,1.915,2.515): ax.text(x0+xx,5.18,"✂",ha="center",va="center",fontsize=9,color=TEXT2)
ax.text(x0+0.35,5.0,"i",ha="right",va="center",fontsize=11,color=CI,fontweight="bold")
ax.text(x0+0.35,4.5,"j",ha="right",va="center",fontsize=11,color=CJ,fontweight="bold")

# ③ biotin fill-in
ax.text(2*CW-0.1,4.8,"›",ha="center",va="center",fontsize=24,color=LINE)
head(0,2,3,"Mark the ends","fill in the cut ends with|biotin-labelled nucleotides")
x0=2*CW; prot(x0+1.9,4.75,0.4,0.62)
seg(x0+1.15,x0+1.85,5.0,CI); seg(x0+1.15,x0+1.85,4.5,CJ)
seg(x0+1.98,x0+2.6,5.0,CI); seg(x0+1.98,x0+2.6,4.5,CJ)
for xx,yy in [(1.85,5.0),(1.85,4.5),(1.98,5.0),(1.98,4.5)]: biotin(x0+xx+(-0.02 if xx<1.9 else 0.02),yy)
ax.text(x0+2.85,4.75,"B = biotin",ha="left",va="center",fontsize=8.5,color=CRIMSON)

# ④ proximity ligation
head(1,0,4,"Proximity ligation","ends close in space are joined|→ a chimeric i–j junction")
x0=0.0; prot(x0+1.9,1.4,0.42,0.75)
seg(x0+1.05,x0+1.9,1.65,CI); seg(x0+1.9,x0+2.75,1.65,CJ)
biotin(x0+1.9,1.65)
seg(x0+1.05,x0+1.9,1.15,CJ); seg(x0+1.9,x0+2.75,1.15,CI); biotin(x0+1.9,1.15)
ax.text(x0+0.9,1.65,"i",ha="right",va="center",fontsize=11,color=CI,fontweight="bold")
ax.text(x0+2.9,1.65,"j",ha="left",va="center",fontsize=11,color=CJ,fontweight="bold")

# ⑤ shear + pull-down + paired-end sequencing
ax.text(CW-0.1,1.5,"›",ha="center",va="center",fontsize=24,color=LINE)
head(1,1,5,"Pull down & sequence","reverse crosslinks, shear, keep|biotin junctions, paired-end reads")
x0=CW
seg(x0+0.65,x0+1.75,1.4,CI); seg(x0+1.75,x0+2.85,1.4,CJ); biotin(x0+1.75,1.4)
ax.add_patch(mp.Circle((x0+1.75,1.8),0.17,fc=TEXT2,ec="none",zorder=4))
ax.plot([x0+1.75,x0+1.75],[1.48,1.64],color=TEXT2,lw=1.2,zorder=3)
ax.text(x0+2.0,1.83,"streptavidin bead",ha="left",va="center",fontsize=8.5,color=TEXT2)
ax.annotate("",(x0+1.3,1.12),(x0+0.68,1.12),arrowprops=dict(arrowstyle="-|>",color=CI,lw=1.6))
ax.annotate("",(x0+2.2,1.12),(x0+2.82,1.12),arrowprops=dict(arrowstyle="-|>",color=CJ,lw=1.6))
ax.text(x0+0.68,0.9,"read 1",ha="left",va="center",fontsize=8.5,color=CI)
ax.text(x0+2.82,0.9,"read 2",ha="right",va="center",fontsize=8.5,color=CJ)

# ⑥ map the pair → +1 in the contact matrix
ax.text(2*CW-0.1,1.5,"›",ha="center",va="center",fontsize=24,color=LINE)
head(1,2,6,"Count contacts","the two reads map to i and j|→ +1 in cell (i, j)")
x0=2*CW; M0x,M0y,S=x0+1.35,0.0,1.75; n=7; c=S/n
rng=np.random.default_rng(2)
for a in range(n):
    for b in range(n):
        d=abs(a-b); v=0.85 if d==0 else (0.45 if d==1 else 0.12+0.06*rng.random())
        ax.add_patch(mp.Rectangle((M0x+b*c,M0y+(n-1-a)*c),c,c,fc=tint(v*0.75),ec="white",lw=0.8,zorder=2))
ia,jb=1,5
for (a,b) in [(ia,jb),(jb,ia)]:
    ax.add_patch(mp.Rectangle((M0x+b*c,M0y+(n-1-a)*c),c,c,fc=CRIMSON,ec="white",lw=0.8,zorder=3))
ax.text(M0x+jb*c+c/2,M0y+(n-1-ia)*c+c/2,"+1",ha="center",va="center",fontsize=8,color="white",fontweight="bold",zorder=4)
ax.text(M0x+ia*c+c/2,M0y+S+0.12,"i",ha="center",va="bottom",fontsize=10,color=CI,fontweight="bold")
ax.text(M0x+jb*c+c/2,M0y+S+0.12,"j",ha="center",va="bottom",fontsize=10,color=CJ,fontweight="bold")
ax.text(M0x-0.1,M0y+(n-1-ia)*c+c/2,"i",ha="right",va="center",fontsize=10,color=CI,fontweight="bold")
ax.text(M0x-0.1,M0y+(n-1-jb)*c+c/2,"j",ha="right",va="center",fontsize=10,color=CJ,fontweight="bold")
ax.text(x0+0.08,1.3,"contact",ha="left",va="center",fontsize=8.8,color=TEXT,fontweight="semibold")
ax.text(x0+0.08,1.05,"matrix",ha="left",va="center",fontsize=8.8,color=TEXT,fontweight="semibold")
ax.text(x0+0.08,0.75,"genome bins",ha="left",va="center",fontsize=8.5,color=TEXT2)
ax.text(x0+0.08,0.52,"× genome bins",ha="left",va="center",fontsize=8.5,color=TEXT2)

ax.plot([0,W],[-0.35,-0.35],color="#E3E6EB",lw=0.9)
ax.text(0.02,-0.48,"Millions of read pairs fill the matrix: each cell counts how often two genome bins were found together. Schematic.",
        ha="left",va="top",fontsize=9.5,color=TEXT)
header(ax,"Hi-C · genome-wide chromosome conformation capture","Loci that touch in 3D are ligated into one molecule",dy=4)
fig.savefig(OUT+"hic_principle.png"); print("ok g12 hic principle")
