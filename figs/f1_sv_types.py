from gallery import *

fig, ax = plt.subplots(figsize=(11.5,5.4))
W,G,H=0.9,0.06,0.42          # segment width, gap, height

def seg(x,y,name,rev=False,w=W,empty=False,hot=False,other=False):
    t=0.16                   # arrow-tip length: shows orientation
    if rev: pts=[(x+w,y),(x+t,y),(x,y+H/2),(x+t,y+H),(x+w,y+H)]
    else:   pts=[(x,y),(x+w-t,y),(x+w,y+H/2),(x+w-t,y+H),(x,y+H)]
    if empty:
        ax.add_patch(mp.Polygon(pts,closed=True,fc="white",ec=CRIMSON,lw=1,ls=(0,(2,2))))
        return x+w+G
    fc=CRIMSON if hot else ("white" if other else FAINT)
    ax.add_patch(mp.Polygon(pts,closed=True,fc=fc,ec=LINE if other else "none",lw=0.9))
    ax.text(x+w/2+(0.05 if rev else -0.05),y+H/2,name,ha="center",va="center",fontsize=11,
            color="white" if hot else TEXT2,fontweight="bold" if hot else "semibold")
    return x+w+G

def row(y,items):
    x=0
    for it in items:
        if it=="|": x+=0.35; continue
        name,kw=(it,{}) if isinstance(it,str) else it
        x=seg(x,y,name,**kw)

oth={"other":True}
rows=[("Reference",["A","B","C","D","|",("W",oth),("X",oth)],("chromosome 1: A B C D  ·  chromosome 2: W X","")),
      ("Deletion",["A",("B",{"empty":True}),"C","D"],("B lost → copy number ↓","CNV")),
      ("Duplication",["A",("B",{"hot":True}),("B",{"hot":True}),"C","D"],("B copied → copy number ↑","CNV")),
      ("Insertion",["A","B",("new",{"w":0.78,"hot":True}),"C","D"],("new sequence, often a mobile element (Alu, L1, SVA)","")),
      ("Inversion",["A",("B",{"rev":True,"hot":True}),"C","D"],("B flipped → copy number unchanged","balanced")),
      ("Translocation",["A","B",("X",{"hot":True,"other":True}),"|",("W",oth),("C",{"hot":True}),("D",{"hot":True})],
       ("ends swapped between chromosomes","balanced"))]

for i,(lab,items,(note,kind)) in enumerate(rows):
    y=-i*0.72
    ax.text(-0.15,y+H/2,lab,ha="right",va="center",fontsize=11.5,color=TEXT,
            fontweight="bold" if i==0 else "semibold")
    row(y,items)
    ax.text(6.55,y+H/2+(0.1 if kind else 0),note,ha="left",va="center",fontsize=10,color=TEXT2)
    if kind: ax.text(6.55,y+H/2-0.14,kind,ha="left",va="center",fontsize=10,color=CRIMSON,fontweight="semibold")
    if i==0: ax.plot([-1.6,10.9],[y-0.15,y-0.15],color="#E3E6EB",lw=0.9)

ax.set_xlim(-1.6,10.9); ax.set_ylim(-5*0.72-0.1,0.6)
ax.axis("off")
header(ax,"Structural variants (≥ 50 bp)","Five shapes — only deletions and duplications change copy number")
fig.savefig(OUT+"sv_types.png"); print("ok f1 sv")
