from style import *
import matplotlib.patches as mp
fig=plt.figure(figsize=(12.6,6.2))
gs=fig.add_gridspec(2,2,height_ratios=[0.85,1.25],hspace=0.50,wspace=0.26)

ax=fig.add_subplot(gs[0,:]); ax.axis("off")
ax.set_title("The 16S rRNA gene: conserved regions let you amplify, variable regions let you classify",
             loc="left",x=0,fontsize=12.5)
L=1540
ax.add_patch(mp.Rectangle((0,0.40),1.0,0.26,fc="#DCE2E7",ec="none"))
V=[(69,99),(137,242),(433,497),(576,682),(822,879),(986,1043),(1117,1173),(1243,1294),(1435,1465)]
for i,(a,b) in enumerate(V):
    ax.add_patch(mp.Rectangle((a/L,0.40),(b-a)/L,0.26,fc=TEAL,ec="none"))
    ax.text((a+b)/2/L,0.30,f"V{i+1}",ha="center",fontsize=8.6,color="#2E5E4E")
ax.text(-0.012,0.53,"16S",ha="right",va="center",fontsize=10.5,fontweight="bold")
ax.annotate("",xy=(682/L,0.78),xytext=(433/L,0.78),
            arrowprops=dict(arrowstyle="<|-|>",color=RED,lw=1.8))
ax.text((433+682)/2/L,0.84,"V3–V4  (a common choice)",ha="center",fontsize=9.6,color=RED,fontweight="bold")
ax.text(0.5,0.10,"⚠️  Different primer pairs amplify different taxa with different efficiency. Two studies using different regions are not directly comparable.",
        ha="center",fontsize=9.6,color="#8A3A33")
ax.set_xlim(-0.02,1.02); ax.set_ylim(0,1)

ax=fig.add_subplot(gs[1,0])
levels=["Phylum","Class","Order","Family","Genus","Species","Strain"]
amp=[1.0,1.0,0.97,0.90,0.72,0.22,0.0]
shot=[1.0,1.0,1.0,0.99,0.96,0.85,0.55]
Y=np.arange(len(levels))[::-1]; h=0.36
ax.barh(Y+h/2,amp,h,color="#C9A227",label="16S amplicon")
ax.barh(Y-h/2,shot,h,color=BLUE,label="shotgun")
ax.set_yticks(Y); ax.set_yticklabels(levels,fontsize=9.5)
ax.set_xlabel("fraction of reads confidently classified"); ax.set_xlim(0,1.05)
ax.set_title("Taxonomic resolution",loc="left",fontsize=11.5)
ax.legend(frameon=False,fontsize=9.2,loc="upper center",bbox_to_anchor=(0.5,-0.20),ncol=2)
ax.axhline(Y[4]+0.5,color="#BBB",ls="--",lw=1.1)
ax.text(1.02,Y[6]+0.05,"16S rarely reaches\nspecies or strain level",ha="right",va="center",fontsize=9,
        color="#8A3A33",linespacing=1.4)

ax=fig.add_subplot(gs[1,1]); ax.axis("off")
ax.set_title("Choosing between them",loc="left",fontsize=11.5)
rows=[("16S amplicon","#C9A227",["cheap — hundreds of samples feasible",
        "works at low biomass (PCR amplifies)","bacteria and archaea only",
        "genus level; no gene functions","primer and copy-number bias"]),
      ("Shotgun metagenomics",BLUE,["species and strain resolution",
        "gene and pathway content (HUMAnN)","viruses, fungi, plasmids too",
        "needs high biomass and low host DNA","~10× the cost per sample"])]
y=0.90
for nm,c,items in rows:
    ax.text(0.0,y,nm,fontsize=11,fontweight="bold",color=c); y-=0.075
    for it in items:
        ax.text(0.02,y,"·  "+it,fontsize=9.4,color="#555",va="top"); y-=0.070
    y-=0.075
ax.set_xlim(0,1); ax.set_ylim(0,1)
fig.savefig(OUT+"amplicon_vs_shotgun.png"); print("ok n3")
