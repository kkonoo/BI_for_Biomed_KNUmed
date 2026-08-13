from style import *
import matplotlib.patches as mp
fig,axes=plt.subplots(1,2,figsize=(11.8,4.6),gridspec_kw={"wspace":0.30,"width_ratios":[1,1.15]})

ax=axes[0]
cats=["rRNA","mRNA\n(polyadenylated)","lncRNA / other\nnon-coding","pre-mRNA\n(intronic)"]
cols=["#C1544B","#4A77A8","#5F9EA0","#C9A227"]
total=[0.85,0.03,0.10,0.02]
polyA=[0.02,0.86,0.09,0.03]
ribo =[0.05,0.55,0.29,0.11]
X=np.arange(3); bl=np.zeros(3)
for i,(c,col) in enumerate(zip(cats,cols)):
    v=np.array([total[i],polyA[i],ribo[i]])
    ax.bar(X,v,bottom=bl,color=col,width=0.6,label=c,edgecolor="white",lw=1.2); bl+=v
ax.set_xticks(X); ax.set_xticklabels(["Total RNA\n(no selection)","poly(A)\nselection","rRNA\ndepletion"],fontsize=9.5)
ax.set_ylabel("fraction of sequenced reads"); ax.set_ylim(0,1)
ax.set_title("Where your reads actually go",loc="left",fontsize=12)
ax.legend(frameon=False,fontsize=8.6,loc="upper center",bbox_to_anchor=(0.5,-0.20),ncol=2,
          handletextpad=0.35,columnspacing=0.9)
ax.annotate("without selection ~85% of your\nsequencing budget is rRNA",
            xy=(0.30,0.45),xytext=(1.35,0.62),fontsize=9,color="#8A3A33",ha="center",
            linespacing=1.4,arrowprops=dict(arrowstyle="-|>",color="#C1544B",lw=1.2,
            connectionstyle="arc3,rad=0.25"))

ax=axes[1]
rows=[("poly(A) selection",
       ["cheapest per usable read","clean mRNA quantification","the default for DE studies"],
       ["misses non-polyadenylated RNA\n(histones, many lncRNA, circRNA)","needs intact RNA — fails on FFPE","3′ bias if RNA is degraded"],BLUE),
      ("rRNA depletion (Ribo-Zero etc.)",
       ["works on degraded / FFPE RNA","keeps non-coding and nascent RNA","better for intron / splicing work"],
       ["more reads needed for same mRNA depth","more expensive","higher intronic background"],TEAL)]
ax.axis("off")
y=0.99
for name,pros,cons,c in rows:
    ax.text(0.0,y,name,fontsize=11.5,fontweight="bold",color=c)
    y-=0.080
    for p in pros:
        ax.text(0.02,y,"✓  "+p,fontsize=9.4,color="#2E8B57",va="top"); y-=0.070
    for cc in cons:
        n=cc.count("\n")
        ax.text(0.02,y,"✗  "+cc,fontsize=9.4,color="#B04A42",va="top",linespacing=1.35); y-=0.070*(1+n*0.85)
    y-=0.030
ax.set_xlim(0,1); ax.set_ylim(0,1)
ax.set_title("Choosing between them",loc="left",fontsize=12)
fig.text(0.5,-0.22,"This choice is made at the bench and cannot be undone in analysis. Decide it from the biological question.",
         ha="center",fontsize=9.5,color="#666")
fig.savefig(OUT+"library_selection.png"); print("ok h1")
