from style import *
fig,ax=plt.subplots(figsize=(9.4,5.4))
pts=[("Bulk RNA-seq",1,3e7,GREY,"one averaged 'cell',\nvery deep",34,-46),
     ("Plate-based\n(Smart-seq2/3)",400,2e6,TEAL,"full-length reads,\nisoform-capable,\nlabour-intensive",42,-58),
     ("Droplet\n(10x Genomics)",50000,25000,BLUE,"3′ or 5′ end only,\nthe current default",42,-46),
     ("Combinatorial\nindexing (sci-)",1_000_000,4000,PURPLE,"whole-organism atlases,\nvery shallow per cell",42,-46)]
for nm,n,depth,c,note,dy1,dy2 in pts:
    ax.scatter([n],[depth],s=360,color=c,alpha=0.26,zorder=3)
    ax.scatter([n],[depth],s=85,color=c,zorder=4)
    ax.annotate(nm,(n,depth),xytext=(0,dy1),textcoords="offset points",ha="center",
                fontsize=10.5,fontweight="bold",color=c,linespacing=1.3)
    ax.annotate(note,(n,depth),xytext=(0,dy2),textcoords="offset points",ha="center",
                fontsize=8.8,color="#666",style="italic",linespacing=1.4)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(0.3,6e6); ax.set_ylim(6e2,4e8)
ax.set_xlabel("cells profiled per experiment"); ax.set_ylabel("reads (or UMIs) per cell")
ax.set_title("The single-cell trade-off: how many cells vs. how well",loc="left")
x=np.logspace(2.2,6.3,80)
ax.plot(x,1.3e9/x,color="#CFCFCF",lw=2.0,ls="--",zorder=1)
ax.text(4e5,1.3e9/4e5*2.4,"one fixed sequencing budget",fontsize=9.2,color="#9A9A9A",rotation=-25)
fig.text(0.5,-0.02,"The three single-cell options sit on roughly one budget line. Choosing where on that line you want to be is the first design question.",
         ha="center",fontsize=9.5,color="#666")
fig.savefig(OUT+"platform_tradeoff.png"); print("ok")
