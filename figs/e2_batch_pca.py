from style import *
rng=np.random.default_rng(12)
G=3000
base=rng.normal(6,2,G)
batch_eff=np.where(rng.random(G)<0.45,rng.normal(0,0.9,G),0)
dis_eff=np.where(rng.random(G)<0.15,rng.normal(0,1.0,G),0)
indiv=rng.normal(0,0.09,(4,G))   # person-to-person programs, so groups have spread

def simulate(batch,case):
    X=(base+np.outer(batch,batch_eff)+np.outer(case,dis_eff)
       +rng.normal(0,1,(len(batch),4))@indiv+rng.normal(0,0.6,(len(batch),G)))
    X=X-X.mean(0)
    U,S,_=np.linalg.svd(X,full_matrices=False)
    pc=U*S; var=S**2/np.sum(S**2)*100
    return pc[:,:2]*np.sign(pc[np.argmax(batch),:2]),var

designs=[("Correctable: each batch = 5 case + 5 control",
          np.repeat([0,1],10),np.tile(np.repeat([1,0],5),2)),
         ("Not analysable: batch 1 = case, batch 2 = control",
          np.repeat([0,1],10),np.repeat([1,0],10))]
notes=["PC1 = batch — bigger than the disease effect.\nBut case and control sit in both batches,\nso the disease shows up on PC2 → correctable",
       "Two clean groups — but colour (disease) and\nshape (batch) coincide perfectly.\nDisease or batch? No method can tell"]

fig,axes=plt.subplots(1,2,figsize=(11.6,4.9),gridspec_kw={"wspace":0.22})
for ax,(title,batch,case),note,nc in zip(axes,designs,notes,["#2E5E4E","#8A3A33"]):
    pc,var=simulate(batch,case)
    for b,mk in [(0,"o"),(1,"^")]:
        for c,col in [(1,RED),(0,BLUE)]:
            m=(batch==b)&(case==c)
            ax.scatter(pc[m,0],pc[m,1],s=70,marker=mk,color=col,edgecolor="white",lw=0.8,alpha=0.9,zorder=3)
    ax.set_xlabel(f"PC1 ({var[0]:.0f}% of variance)"); ax.set_ylabel(f"PC2 ({var[1]:.0f}%)")
    ax.set_xticks([]); ax.set_yticks([])
    lim=np.abs(pc).max()*1.25
    ax.set_xlim(-lim,lim); ax.set_ylim(-lim*0.75,lim*0.75)
    ax.set_title(title,loc="left",fontsize=11.5,color=nc)
    ax.text(0.5,-0.14,note,transform=ax.transAxes,ha="center",va="top",fontsize=9.6,color=nc,linespacing=1.45)

from matplotlib.lines import Line2D
h=[Line2D([],[],marker="o",ls="",color=RED,ms=8,label="case"),
   Line2D([],[],marker="o",ls="",color=BLUE,ms=8,label="control"),
   Line2D([],[],marker="o",ls="",color="#777",ms=8,label="batch 1"),
   Line2D([],[],marker="^",ls="",color="#777",ms=8,label="batch 2")]
fig.legend(handles=h,loc="upper center",ncol=4,frameon=False,fontsize=10,bbox_to_anchor=(0.5,1.06))
fig.savefig(OUT+"batch_pca.png"); print("ok e2")
