from style import *
rng=np.random.default_rng(3)
fig,axes=plt.subplots(1,4,figsize=(14.2,4.2),gridspec_kw={"wspace":0.30,"width_ratios":[0.72,1,1,1]})

ax=axes[0]
ax.bar([0,1],[100,200],color=[GREY,RED],width=0.55,edgecolor="white")
ax.set_xticks([0,1]); ax.set_xticklabels(["Ctrl","Treat"])
ax.set_ylabel("bulk expression of GENE-X"); ax.set_ylim(0,260)
ax.text(0.5,0.93,"2.0×\np < 0.001",transform=ax.transAxes,ha="center",va="top",
        fontsize=11,fontweight="bold",color=RED,linespacing=1.4)
ax.set_title("Bulk result",loc="left",fontsize=12)

titles=["① every cell shifts a little","② a subpopulation shifts a lot","③ the mixture changed"]
notes=["a general, modest response","a rare responsive population\n— the interesting case",
       "no cell changed at all\n(see W4 §7.2)"]
cols=[BLUE,TEAL,RED]
for k,(ax,t,note,c) in enumerate(zip(axes[1:],titles,notes,cols)):
    n=600
    if k==0:
        ctrl=rng.lognormal(4.3,0.55,n); trt=rng.lognormal(4.3+np.log(2),0.55,n)
    elif k==1:
        ctrl=rng.lognormal(4.3,0.55,n)
        trt=np.concatenate([rng.lognormal(4.3,0.55,int(n*0.85)),
                            rng.lognormal(4.3+np.log(8),0.45,int(n*0.15))])
    else:
        lo=rng.lognormal(3.6,0.45,n); hi=rng.lognormal(5.3,0.45,n)
        ctrl=np.concatenate([lo[:480],hi[:120]])
        trt =np.concatenate([lo[:170],hi[:430]])
    bins=np.linspace(0,600,60)
    ax.hist(ctrl,bins=bins,color=GREY,alpha=0.72,lw=0,label="Ctrl")
    ax.hist(trt,bins=bins,color=c,alpha=0.55,lw=0,label="Treat")
    ax.set_xlabel("expression per cell"); ax.set_yticks([]); ax.set_xlim(0,600)
    ax.set_title(t,loc="left",fontsize=11.5,color=c)
    ax.legend(frameon=False,fontsize=9,loc="upper right")
    ax.text(0.5,-0.30,note,transform=ax.transAxes,ha="center",va="top",fontsize=9.4,
            color="#555",linespacing=1.4)
    if k==0: ax.set_ylabel("number of cells")
fig.suptitle("The same bulk result. Three completely different biologies.",
             x=0.075,ha="left",fontsize=13.5,fontweight="bold",y=1.02)
fig.savefig(OUT+"bulk_vs_singlecell.png"); print("ok k1")
