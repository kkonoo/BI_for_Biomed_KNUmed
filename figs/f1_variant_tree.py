from gallery import *

fig, ax = plt.subplots(figsize=(12,5.0))
ax.set_xlim(0,1); ax.set_ylim(0.17,0.9); ax.axis("off")
header(ax,"Classes of human genetic variation","Split by size — except tandem repeats, defined by sequence")
G="#9AA1AB"   # connector grey

# root: the one dark box
card(ax,0.5,0.825,0.40,0.085,fc=CRIMSON,ec=CRIMSON)
ax.text(0.5,0.838,"Human genetic variation",ha="center",va="center",color="white",fontsize=13.5,fontweight="bold",zorder=3)
ax.text(0.5,0.806,"differences from the reference genome",ha="center",va="center",color=tint(0.18),fontsize=9.5,zorder=3)

cols=[(0.17,"1","Small variants","< 50 bp  ·  ~4–5 million per genome",
       [("SNV","one base substituted"),("Indel","1–49 bp inserted or deleted")],
       "arrays · short reads"),
      (0.5,"2","Structural variants","≥ 50 bp  ·  ~20–25 thousand per genome",
       [("Deletion · Duplication","copy number changes = CNV"),
        ("Insertion","incl. mobile elements: Alu, L1, SVA"),
        ("Inversion · Translocation","balanced: copy number unchanged")],
       "long reads best · short reads partly"),
      (0.83,"3","Tandem repeats","one unit repeated side by side",
       [("STR (microsatellite)","unit 1–6 bp, e.g. (CAG)n"),("VNTR (minisatellite)","unit ≥ 7 bp")],
       "dedicated callers · long reads")]
W=0.315
ax.plot([0.5,0.5],[0.782,0.745],color=G,lw=1,zorder=1)
ax.plot([0.17,0.83],[0.745,0.745],color=G,lw=1,zorder=1)
for x,num,title,sub,leaves,tech in cols:
    ax.plot([x,x],[0.745,0.71],color=G,lw=1,zorder=1)
    card(ax,x,0.665,W,0.09,fc=tint(0.16),ec=tint(0.16))
    ax.scatter([x-W/2+0.03],[0.665],s=420,color=CRIMSON,zorder=3)
    ax.text(x-W/2+0.03,0.665,num,ha="center",va="center",color="white",fontsize=10,fontweight="bold",zorder=4)
    ax.text(x-W/2+0.06,0.68,title,ha="left",va="center",color=CRIMSON,fontsize=12.5,fontweight="bold",zorder=3)
    ax.text(x-W/2+0.06,0.646,sub,ha="left",va="center",color=TEXT2,fontsize=8.9,zorder=3)
    lx=x-W/2+0.02; ys=[0.555-i*0.098 for i in range(len(leaves))]
    ax.plot([lx,lx],[0.62,ys[-1]],color=G,lw=1,zorder=1)
    for y,(name,desc) in zip(ys,leaves):
        ax.plot([lx,lx+0.025],[y,y],color=G,lw=1,zorder=1)
        card(ax,lx+0.025+0.13,y,0.26,0.078)
        ax.text(lx+0.04,y+0.014,name,ha="left",va="center",fontsize=10.5,fontweight="semibold",color=TEXT,zorder=3)
        ax.text(lx+0.04,y-0.017,desc,ha="left",va="center",fontsize=9,color=TEXT2,zorder=3)
    ax.text(x-W/2,0.235,"detected by",ha="left",va="center",fontsize=8.8,color=TEXT2)
    ax.text(x-W/2,0.2,tech,ha="left",va="center",fontsize=10,color=TEXT,fontweight="semibold")
ax.plot([0.015,0.985],[0.27,0.27],color="#E3E6EB",lw=0.9)
fig.savefig(OUT+"variant_tree.png"); print("ok f1 tree")
