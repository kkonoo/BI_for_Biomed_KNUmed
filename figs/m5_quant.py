from style import *
import matplotlib.patches as mp
fig,axes=plt.subplots(1,2,figsize=(12.8,4.9),gridspec_kw={"wspace":0.26,"width_ratios":[1.3,1]})

ax=axes[0]; ax.axis("off")
ax.set_title("When are the samples combined?",loc="left",fontsize=12.5)
X=[0.35,0.545,0.74,0.935]
for x,s in zip(X,["cells","lysate","peptides","LC-MS run"]):
    ax.text(x,0.99,s,ha="center",fontsize=9.4,color="#777")
rows=[("SILAC",0,"#2E8B57","metabolic labelling — earliest mixing,\nso every later step is shared"),
      ("TMT / iTRAQ",2,"#C9A227","chemical labelling after digestion;\nup to 18 samples in one run"),
      ("Label-free (LFQ)",3,BLUE,"never mixed — every sample is its own\nrun, so every step can drift")]
y=0.80
for nm,mixat,c,note in rows:
    ax.text(0.0,y,nm,fontsize=11,fontweight="bold",color=c,va="center")
    ax.text(0.0,y-0.095,note,fontsize=9.1,color="#666",va="top",linespacing=1.45)
    for j,x in enumerate(X):
        if j<mixat:
            ax.add_patch(mp.Circle((x-0.033,y),0.018,fc=c,ec="none"))
            ax.add_patch(mp.Circle((x+0.033,y),0.018,fc="#9AA5B1",ec="none"))
        else:
            ax.add_patch(mp.Ellipse((x,y),0.086,0.046,fc="#EDEFF1",ec="#BFC7CE",lw=1.1))
            ax.add_patch(mp.Circle((x-0.016,y),0.016,fc=c,ec="none"))
            ax.add_patch(mp.Circle((x+0.016,y),0.016,fc="#9AA5B1",ec="none"))
        if j<3: ax.annotate("",xy=(x+0.078,y),xytext=(x+0.058,y),
                            arrowprops=dict(arrowstyle="-|>",color="#D5D5D5",lw=1.3))
    xm=X[mixat]
    ax.annotate("",xy=(xm,y+0.048),xytext=(xm,y+0.115),
                arrowprops=dict(arrowstyle="-|>",color=c,lw=1.7))
    ax.text(xm,y+0.125,"mix here",ha="center",fontsize=8.8,color=c,fontweight="bold")
    y-=0.315
ax.text(0.0,-0.10,"→  The earlier you mix, the less technical variation accumulates.",
        fontsize=9.8,color="#2E5E4E",fontweight="bold",va="top")
ax.set_xlim(-0.02,1.0); ax.set_ylim(-0.18,1.06)

ax=axes[1]
true=np.logspace(0,np.log10(12),80)
ax.plot(true,true,color="#999",ls="--",lw=1.6,label="no compression (ideal)")
ax.plot(true,1+(true-1)*0.86,color=TEAL,lw=2.3,label="SPS-MS3 (largely corrected)")
ax.plot(true,1+(true-1)/(1+0.34*(true-1)),color=RED,lw=2.5,label="MS2-based TMT (observed)")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("true abundance ratio"); ax.set_ylabel("measured ratio")
ax.set_xticks([1,2,4,8,12]); ax.set_xticklabels(["1","2","4","8","12"])
ax.set_yticks([1,2,4,8,12]); ax.set_yticklabels(["1","2","4","8","12"])
ax.set_title("⚠️ TMT ratio compression",loc="left",fontsize=12.5,color=RED)
ax.legend(frameon=False,fontsize=9.2,loc="upper left")
ax.text(0.985,0.06,"Peptides co-isolated with your target contribute\nto the same reporter ions, diluting the signal.\nA real 8× change can be reported as 3×.",
        transform=ax.transAxes,ha="right",fontsize=9.2,color="#333",linespacing=1.5,
        bbox=dict(boxstyle="round,pad=0.42",fc="#FDF3F2",ec="#E8C4C0"))
fig.savefig(OUT+"quantification.png"); print("ok")
