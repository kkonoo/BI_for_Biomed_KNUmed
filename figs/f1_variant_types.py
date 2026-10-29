from gallery import *

fig, ax = plt.subplots(figsize=(10,4.3))
ref=list("ACGTTGCATGCAGTTACG")
rows=[
 ("Reference", ref, None),
 ("SNV", list("ACGTTGCATGCTGTTACG"), [11]),
 ("Insertion", list("ACGTTGCATGCAG")+list("TTA")+list("TTACG"), list(range(13,16))),
 ("Deletion", list("ACGTTGCA")+["-"]*3+list("AGTTACG"), [8,9,10]),
]
notes=["","1 base substituted","3 bases inserted","3 bases deleted"]
S,B=0.055,0.048              # step and box width
for i,(label,seq,hl) in enumerate(rows):
    y=len(rows)-i-1
    for j,b in enumerate(seq):
        on=hl is not None and j in hl
        if b=="-":
            ax.add_patch(mp.Rectangle((j*S,y),B,0.5,fc="white",ec=CRIMSON,lw=1,ls=(0,(2,2))))
            continue
        ax.add_patch(mp.Rectangle((j*S,y),B,0.5,fc=CRIMSON if on else FAINT,ec="none"))
        ax.text(j*S+B/2,y+0.25,b,ha="center",va="center",fontsize=11,
                color="white" if on else TEXT2,fontweight="bold" if on else "normal")
    ax.text(-0.03,y+0.25,label,ha="right",va="center",fontsize=11.5,
            fontweight="bold" if i==0 else "semibold",color=TEXT)
    if notes[i]: ax.text(len(seq)*S+0.03,y+0.25,notes[i],ha="left",va="center",fontsize=10,color=TEXT2)

ax.set_xlim(-0.32,1.45); ax.set_ylim(-0.15,len(rows)-0.4)
ax.axis("off")
header(ax,"Small variants (< 50 bp)","SNVs and indels change a few bases against the reference")
fig.savefig(OUT+"small_variants.png"); print("ok f1")
