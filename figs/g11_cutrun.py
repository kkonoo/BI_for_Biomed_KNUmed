from gallery import *

# CUT&RUN and CUT&Tag: an antibody-tethered enzyme cuts only next to the target, in native nuclei
W,Hh=11.6,6.6
fig,ax=plt.subplots(figsize=(11.5,11.5*Hh/W+0.9))
ax.set_xlim(0,W); ax.set_ylim(-0.85,Hh-1.0); ax.set_aspect("equal"); ax.axis("off")
NUC="#D5D9DF"; DARK="#5E6670"; L0=1.75; CW=(W-L0)/4

def dna(x0,x1,y,lw=1.5,c=TEXT2,z=1): ax.plot([x0,x1],[y,y],color=c,lw=lw,solid_capstyle="butt",zorder=z)
def nuc(x,y,mark=False,r=0.12):
    ax.add_patch(mp.Circle((x,y),r,fc=NUC,ec="#A6ADB7",lw=0.8,zorder=2))
    if mark:
        ax.plot([x,x],[y+r,y+r+0.07],color=CRIMSON,lw=1.2,zorder=3)
        ax.add_patch(mp.Circle((x,y+r+0.11),0.05,fc=CRIMSON,ec="none",zorder=3))
def antibody(x,y):                     # arms down onto the mark at (x,y)
    for a,b in [((x,y+0.12),(x,y+0.3)),((x,y+0.12),(x-0.09,y)),((x,y+0.12),(x+0.09,y))]:
        ax.plot([a[0],b[0]],[a[1],b[1]],color=DARK,lw=2,solid_capstyle="round",zorder=4)
def mnase(x,y):
    ax.add_patch(mp.Ellipse((x,y),0.6,0.27,fc=tint(0.45),ec="none",zorder=5))
    ax.text(x,y,"MNase",ha="center",va="center",fontsize=7.5,color=TEXT,zorder=6)
def tn5(x,y,s=0.8):
    for dx in (-0.075*s,0.075*s): ax.add_patch(mp.Circle((x+dx,y),0.085*s,fc=CRIMSON,ec="none",zorder=5))
    ax.plot([x-0.12*s,x-0.24*s],[y+0.05*s,y+0.14*s],color=tint(0.55),lw=2.4,zorder=4)
    ax.plot([x+0.12*s,x+0.24*s],[y+0.05*s,y+0.14*s],color=tint(0.55),lw=2.4,zorder=4)
def adapter(x,y,d):
    ax.add_patch(mp.Rectangle((x if d>0 else x-0.15,y-0.045),0.15,0.09,fc=tint(0.55),ec="none",zorder=4))
def chrom(x0,y,cuts=None,adapters=False):
    dna(x0+0.1,x0+2.25,y)
    for xx in (0.35,0.65,1.55,1.85): nuc(x0+xx,y)
    nuc(x0+1.1,y,mark=True)
    if cuts:
        for xx in (0.86,1.34):
            ax.plot([x0+xx,x0+xx],[y-0.07,y+0.07],color="white",lw=4,zorder=3)
        if adapters:
            for xx,d in [(0.86,-1),(0.86,1),(1.34,-1),(1.34,1)]: adapter(x0+xx+0.02*d,y,d)
def cap(col,y,*lines):
    x0=L0+col*CW
    for k,l in enumerate(lines): ax.text(x0+0.08,y-k*0.25,l,ha="left",va="top",fontsize=8.8,color=TEXT2)

# column heads (shared)
heads=["Antibody binds","Enzyme docks","Activate & cut","Collect & sequence"]
for i,t in enumerate(heads):
    x0=L0+i*CW
    ax.add_patch(mp.Circle((x0+0.2,4.75),0.15,fc=CRIMSON,ec="none"))
    ax.text(x0+0.2,4.75,str(i+1),ha="center",va="center",fontsize=10,fontweight="bold",color="white")
    ax.text(x0+0.43,4.75,t,ha="left",va="center",fontsize=11,fontweight="semibold",color=TEXT)
    if i: ax.text(x0-0.08,3.55,"›",ha="center",va="center",fontsize=22,color=LINE)
    if i: ax.text(x0-0.08,1.25,"›",ha="center",va="center",fontsize=22,color=LINE)
ax.text(0.05,5.25,"native cells or nuclei, permeabilized and bound to beads · no crosslinking, no sonication",
        ha="left",va="center",fontsize=9,color=TEXT2,style="italic")

# row labels
for y,t,sub in [(3.6,"CUT&RUN","MNase"),(1.3,"CUT&Tag","Tn5")]:
    ax.text(0.05,y+0.12,t,ha="left",va="center",fontsize=12,fontweight="bold",color=TEXT)
    ax.text(0.05,y-0.16,"enzyme: "+sub,ha="left",va="center",fontsize=9,color=TEXT2)
ax.plot([0.05,W],[2.4,2.4],color="#E3E6EB",lw=0.9)

# --- CUT&RUN row (y=3.4)
y=3.4
x0=L0;        chrom(x0,y); antibody(x0+1.1,y+0.29)
cap(0,y-0.3,"antibody binds the","target inside the nucleus")
x0=L0+CW;     chrom(x0,y); antibody(x0+1.1,y+0.29); mnase(x0+1.1,y+0.72)
ax.text(x0+1.48,y+0.72,"pA-MNase",ha="left",va="center",fontsize=8.5,color=TEXT2)
cap(1,y-0.3,"protein A/G–MNase","binds the antibody")
x0=L0+2*CW;   chrom(x0,y,cuts=True); antibody(x0+1.1,y+0.29); mnase(x0+1.1,y+0.72)
ax.text(x0+0.25,y+0.6,"Ca²⁺",ha="center",va="center",fontsize=10,color=CRIMSON,fontweight="bold")
cap(2,y-0.3,"Ca²⁺ switches MNase on","→ cuts on both sides")
x0=L0+3*CW
ax.add_patch(mp.Arc((x0+0.3,y+0.25),1.3,1.2,theta1=-50,theta2=50,color=LINE,lw=1.2,ls=(0,(3,2))))
ax.text(x0+0.55,y+0.95,"nucleus",ha="center",va="center",fontsize=8,color=TEXT2,style="italic")
dna(x0+1.3,x0+1.78,y+0.25,lw=1.6); nuc(x0+1.54,y+0.25,mark=True)
ax.annotate("",(x0+1.25,y+0.25),(x0+0.75,y+0.05),arrowprops=dict(arrowstyle="-|>",color=TEXT2,lw=1.1))
cap(3,y-0.3,"the cut fragment diffuses","out → library → sequence")

# --- CUT&Tag row (y=1.1)
y=1.1
x0=L0;        chrom(x0,y); antibody(x0+1.1,y+0.29)
cap(0,y-0.3,"antibody binds the","target inside the nucleus")
x0=L0+CW;     chrom(x0,y); antibody(x0+1.1,y+0.29); tn5(x0+1.1,y+0.72)
ax.text(x0+1.42,y+0.72,"pA-Tn5",ha="left",va="center",fontsize=8.5,color=TEXT2)
cap(1,y-0.3,"protein A–Tn5, loaded","with adapters, binds")
x0=L0+2*CW;   chrom(x0,y,cuts=True,adapters=True); antibody(x0+1.1,y+0.29); tn5(x0+1.1,y+0.72)
ax.text(x0+0.25,y+0.6,"Mg²⁺",ha="center",va="center",fontsize=10,color=CRIMSON,fontweight="bold")
cap(2,y-0.3,"Mg²⁺ → cut and paste","adapters next to the target")
x0=L0+3*CW
dna(x0+0.35,x0+1.55,y+0.25,lw=2); adapter(x0+0.35,y+0.25,-1); adapter(x0+1.55,y+0.25,1)
ax.annotate("",(x0+0.95,y+0.05),(x0+0.4,y+0.05),arrowprops=dict(arrowstyle="-|>",color=CRIMSON,lw=1.4))
ax.annotate("",(x0+1.0,y+0.05),(x0+1.55,y+0.05),arrowprops=dict(arrowstyle="-|>",color=CRIMSON,lw=1.4))
ax.text(x0+1.85,y+0.25,"PCR",ha="left",va="center",fontsize=8.5,color=TEXT2)
cap(3,y-0.3,"fragments already carry","adapters → PCR → sequence")

ax.plot([0.05,W],[0.05,0.05],color="#E3E6EB",lw=0.9)
ax.text(0.05,-0.12,"Only the target region is cut or tagged → little background → thousands instead of millions of cells, "
        "and a fraction of ChIP's reads.",ha="left",va="top",fontsize=9.5,color=TEXT,fontweight="semibold")
ax.text(0.05,-0.47,"Caveat: CUT&Tag uses Tn5, which prefers open chromatin → some bias toward accessible regions. Schematic.",
        ha="left",va="top",fontsize=9.2,color=TEXT2)
header(ax,"CUT&RUN and CUT&Tag","Cut only at the target, inside the nucleus",dy=4)
fig.savefig(OUT+"cutrun_principle.png"); print("ok g11 cutrun")
