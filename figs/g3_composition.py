from style import *
rng=np.random.default_rng(13)
fig,axes=plt.subplots(1,3,figsize=(13.2,4.4),gridspec_kw={"wspace":0.34,"width_ratios":[1,1,1.05]})

ct=["T cell","B cell","Monocyte","NK","Neutrophil"]
cols=["#4A77A8","#5F9EA0","#C9A227","#9B8EC4","#C1544B"]
ctrl=np.array([0.22,0.09,0.10,0.06,0.53])
case=np.array([0.14,0.06,0.13,0.04,0.63])
ax=axes[0]
bl=np.zeros(2)
for i,(c,col) in enumerate(zip(ct,cols)):
    v=np.array([ctrl[i],case[i]])
    ax.bar(["Controls","Cases"],v,bottom=bl,color=col,width=0.55,label=c,edgecolor="white",lw=1.2)
    bl+=v
ax.set_ylabel("cell type proportion"); ax.set_ylim(0,1)
ax.set_title("① Cases have different\n    blood cell composition",loc="left",fontsize=11.5)
ax.legend(frameon=False,fontsize=8.5,loc="upper center",bbox_to_anchor=(0.5,-0.10),ncol=3,
          handletextpad=0.35,columnspacing=0.8)

ax=axes[1]
x=np.arange(5)
meth=np.array([0.22,0.31,0.78,0.28,0.83])
ax.bar(x,meth,color=cols,width=0.62,edgecolor="white",lw=1.2)
ax.set_xticks(x); ax.set_xticklabels(ct,rotation=30,ha="right",fontsize=9)
ax.set_ylabel("β value at CpG cg12345678"); ax.set_ylim(0,1)
ax.set_title("② This CpG is methylated\n    only in some cell types",loc="left",fontsize=11.5)

ax=axes[2]
sd=0.035
c1=np.dot(ctrl,meth)+rng.normal(0,sd,140)
c2=np.dot(case,meth)+rng.normal(0,sd,140)
parts=ax.violinplot([c1,c2],positions=[0,1],widths=0.62,showmeans=True)
for b,col in zip(parts["bodies"],[GREY,RED]):
    b.set_facecolor(col); b.set_alpha(0.55); b.set_edgecolor("none")
for k in ["cmins","cmaxes","cbars","cmeans"]: parts[k].set_color("#666")
ax.set_xticks([0,1]); ax.set_xticklabels(["Controls","Cases"])
ax.set_ylabel("β value")
ax.set_title("③ Apparent differential\n    methylation — p = 3×10⁻⁹",loc="left",fontsize=11.5,color=RED)
ax.text(0.5,0.06,"No cell has changed its methylation.\nOnly the mixture changed.",
        transform=ax.transAxes,ha="center",fontsize=10,color="#8A3A33",fontweight="bold",
        linespacing=1.5,bbox=dict(boxstyle="round,pad=0.45",fc="#FDF3F2",ec="#E8C4C0"))
fig.suptitle("⚠️  Cell-type composition is the dominant confounder in tissue methylation studies",
             x=0.075,ha="left",fontsize=13,fontweight="bold",y=1.04)
fig.savefig(OUT+"cell_composition_confounder.png"); print("ok g3")
