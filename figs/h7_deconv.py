from style import *
import matplotlib.patches as mp
rng=np.random.default_rng(27)
fig,axes=plt.subplots(1,3,figsize=(13.2,4.4),gridspec_kw={"wspace":0.34})
ct=["Neuron","Astrocyte","Microglia","Oligodendrocyte"]
cols=[BLUE,TEAL,"#C9A227",PURPLE]

ax=axes[0]
prof=np.array([[9.2,1.1,0.4,1.9],[1.0,8.6,0.6,2.2],[0.3,0.8,9.0,0.5],[1.5,2.0,0.7,8.8]])
im=ax.imshow(prof,cmap="Blues",aspect="auto")
ax.set_xticks(range(4)); ax.set_xticklabels(["gene\n1","gene\n2","gene\n3","gene\n4"],fontsize=9)
ax.set_yticks(range(4)); ax.set_yticklabels(ct,fontsize=9.5)
for i in range(4):
    for j in range(4):
        ax.text(j,i,f"{prof[i,j]:.1f}",ha="center",va="center",fontsize=9,
                color="white" if prof[i,j]>5 else "#333")
ax.set_title("① Each cell type has its\n    own expression profile",loc="left",fontsize=11.5)

ax=axes[1]
p_ctrl=np.array([0.42,0.26,0.07,0.25])
p_case=np.array([0.31,0.28,0.16,0.25])
bl=np.zeros(2)
for i,(c,col) in enumerate(zip(ct,cols)):
    v=np.array([p_ctrl[i],p_case[i]])
    ax.bar(["Controls","Cases"],v,bottom=bl,color=col,width=0.55,label=c,edgecolor="white",lw=1.2)
    bl+=v
ax.set_ylabel("cell type proportion"); ax.set_ylim(0,1)
ax.set_title("② Disease changes the mixture\n    (neuron loss, gliosis)",loc="left",fontsize=11.5)
ax.legend(frameon=False,fontsize=8.6,loc="upper center",bbox_to_anchor=(0.5,-0.09),ncol=2,
          handletextpad=0.35,columnspacing=0.9)

ax=axes[2]
genes=["gene 1\n(neuronal)","gene 2","gene 3\n(microglial)","gene 4"]
b_ctrl=p_ctrl@prof; b_case=p_case@prof
X=np.arange(4); w=0.36
ax.bar(X-w/2,b_ctrl,w,color=GREY,label="Controls",edgecolor="white")
ax.bar(X+w/2,b_case,w,color=RED,alpha=0.78,label="Cases",edgecolor="white")
ax.set_xticks(X); ax.set_xticklabels(genes,fontsize=8.8)
ax.set_ylabel("bulk expression")
ax.set_title("③ ⚠️ Bulk result",loc="left",fontsize=11.5,color=RED)
ax.legend(frameon=False,fontsize=9,loc="upper right")
for i,(a,b) in enumerate(zip(b_ctrl,b_case)):
    fc=b/a
    ax.text(i,max(a,b)+0.16,f"{fc:.2f}×",ha="center",fontsize=9,fontweight="bold",
            color=RED if abs(np.log2(fc))>0.25 else "#888")
ax.set_ylim(0,max(b_ctrl.max(),b_case.max())*1.32)
ax.text(0.5,0.02,"'Neuronal genes down, microglial genes up'\n— and not one cell changed its expression.",
        transform=ax.transAxes,ha="center",fontsize=9.2,color="#8A3A33",fontweight="bold",linespacing=1.5,
        bbox=dict(boxstyle="round,pad=0.4",fc="#FDF3F2",ec="#E8C4C0"))
fig.suptitle("⚠️  Bulk RNA-seq cannot distinguish a change in expression from a change in cell composition",
             x=0.075,ha="left",fontsize=12.5,fontweight="bold",y=1.05)
fig.savefig(OUT+"deconvolution.png"); print("ok h7")
