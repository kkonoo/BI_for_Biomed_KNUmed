from style import *
import matplotlib.patches as mp

fig, ax = plt.subplots(figsize=(10,5.2))
BASE={'A':"#4A77A8",'C':"#C1544B",'G':"#B5892A",'T':"#5F9EA0",'-':LGREY}
ref=list("ACGTTGCATGCAGTTACG")
rows=[
 ("Reference", ref, None),
 ("SNV\nsingle nucleotide variant", list("ACGTTGCATGCTGTTACG"), [11]),
 ("Insertion", list("ACGTTGCATGCAG")+list("TTA")+list("TTACG"), list(range(13,16))),
 ("Deletion", list("ACGTTGCA")+['-']*3+list("GTTACG"), [8,9,10]),
]
w,h=0.052,0.5
for i,(label,seq,hl) in enumerate(rows):
    y=len(rows)-i-1
    for j,b in enumerate(seq):
        fc=BASE.get(b,LGREY)
        alpha=1.0 if (hl and j in hl) else 0.42
        ax.add_patch(mp.FancyBboxPatch((j*0.055,y),0.05,h,boxstyle="round,pad=0.002,rounding_size=0.01",
                     fc=fc,ec="white",lw=1.2,alpha=alpha))
        ax.text(j*0.055+0.025,y+h/2,b,ha="center",va="center",color="white",
                fontsize=10,fontweight="bold" if (hl and j in hl) else "normal")
    ax.text(-0.03,y+h/2,label,ha="right",va="center",fontsize=10.5)
    if hl:
        ax.annotate("",xy=(hl[0]*0.055+ (len(hl)*0.055)/2 -0.0025, y+h+0.06),
                    xytext=(hl[0]*0.055+(len(hl)*0.055)/2-0.0025, y+h+0.22),
                    arrowprops=dict(arrowstyle="-|>",color=NAVY,lw=1.4))

# SV / CNV panel
yb=-1.35
ax.text(-0.03,yb+0.55,"Structural variant\n> 50 bp",ha="right",va="center",fontsize=10.5)
segs=[(0.0,0.30,TEAL,"A"),(0.30,0.22,PURPLE,"B"),(0.52,0.30,OLIVE,"C")]
for x0,ww,c,lab in segs:
    ax.add_patch(mp.Rectangle((x0,yb+0.62),ww,0.26,fc=c,ec="white",lw=1.5))
    ax.text(x0+ww/2,yb+0.75,lab,ha="center",va="center",color="white",fontweight="bold")
ax.text(0.86,yb+0.75,"reference",ha="left",va="center",fontsize=9,color="#666")
inv=[(0.0,0.30,TEAL,"A"),(0.30,0.22,PURPLE,"B"),(0.52,0.22,PURPLE,"B"),(0.74,0.30,OLIVE,"C")]
for x0,ww,c,lab in inv:
    ax.add_patch(mp.Rectangle((x0,yb+0.18),ww,0.26,fc=c,ec="white",lw=1.5,
                 hatch="///" if lab=="B" else None))
    ax.text(x0+ww/2,yb+0.31,lab,ha="center",va="center",color="white",fontweight="bold")
ax.text(1.06,yb+0.31,"duplication (CNV)",ha="left",va="center",fontsize=9,color=NAVY)

ax.set_xlim(-0.42,1.42); ax.set_ylim(yb-0.05,len(rows)+0.15)
ax.axis("off")
ax.set_title("Classes of human genetic variation",loc="left",x=-0.1)
fig.text(0.5,-0.02,"Highlighted bases mark the difference from the reference. Hatched block = a duplicated segment (copy number variant).",
         ha="center",fontsize=9,color="#666")
fig.savefig(OUT+"variant_types.png")
print("ok")
