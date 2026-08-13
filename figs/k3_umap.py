from style import *
rng=np.random.default_rng(77)
def make(n,cx,cy,sx,sy,rot=0):
    p=np.stack([rng.normal(0,sx,n),rng.normal(0,sy,n)],1)
    c,s=np.cos(rot),np.sin(rot)
    p=p@np.array([[c,-s],[s,c]])
    return p+np.array([cx,cy])

pops=[("CD4 T",1400,(-3.2,3.0),(1.35,0.85),0.5,BLUE),
      ("CD8 T",900,(-1.0,4.6),(0.95,0.70),0.3,"#7FA8D0"),
      ("NK",420,(1.4,5.2),(0.62,0.55),0,"#9B8EC4"),
      ("B",760,(3.6,1.6),(0.85,0.72),0,TEAL),
      ("Monocyte",1150,(2.4,-2.6),(1.15,0.95),-0.4,"#C9A227"),
      ("DC",180,(4.9,-3.9),(0.45,0.40),0,"#E08A3C"),
      ("Platelet",130,(-4.6,-3.2),(0.42,0.38),0,"#C1544B")]
X=[];lab=[]
for nm,n,(cx,cy),(sx,sy),r,c in pops:
    X.append(make(n,cx,cy,sx,sy,r)); lab += [nm]*n
X=np.vstack(X)

fig,axes=plt.subplots(1,3,figsize=(14.0,4.7),gridspec_kw={"wspace":0.14})
ax=axes[0]
i0=0
for nm,n,_,_,_,c in pops:
    ax.scatter(X[i0:i0+n,0],X[i0:i0+n,1],s=4,color=c,lw=0,alpha=0.72)
    cx,cy=X[i0:i0+n].mean(0)
    ax.text(cx,cy,nm,fontsize=9.5,fontweight="bold",ha="center",va="center",
            color="#222",bbox=dict(boxstyle="round,pad=0.22",fc="white",ec="none",alpha=0.78))
    i0+=n
ax.set_title("Annotated UMAP",loc="left",fontsize=12)
for a in axes: a.set_xticks([]); a.set_yticks([]); a.set_xlabel("UMAP 1"); a.set_ylabel("UMAP 2")

# resolution sweep
from itertools import product
def kmeans(X,k,seed):
    r=np.random.default_rng(seed)
    C=X[r.choice(len(X),k,replace=False)]
    for _ in range(28):
        d=((X[:,None,:]-C[None])**2).sum(2); a=d.argmin(1)
        for j in range(k):
            if (a==j).any(): C[j]=X[a==j].mean(0)
    return a
pal=plt.get_cmap("tab20")
for ax,k,res in zip(axes[1:],[4,14],[0.1,1.6]):
    a=kmeans(X.copy(),k,5)
    ax.scatter(X[:,0],X[:,1],s=4,c=[pal(j%20) for j in a],lw=0,alpha=0.75)
    ax.set_title(f"resolution = {res}  →  {k} clusters",loc="left",fontsize=12)
axes[1].text(0.5,0.02,"CD4 and CD8 T are indistinguishable here",transform=axes[1].transAxes,
             ha="center",fontsize=9.3,color="#8A3A33",
             bbox=dict(boxstyle="round,pad=0.35",fc="#FDF3F2",ec="#E8C4C0"))
axes[2].text(0.5,0.02,"monocytes split into four —\nis that biology or noise?",transform=axes[2].transAxes,
             ha="center",fontsize=9.3,color="#8A3A33",linespacing=1.4,
             bbox=dict(boxstyle="round,pad=0.35",fc="#FDF3F2",ec="#E8C4C0"))
fig.suptitle("Same data, same algorithm — one arbitrary parameter",
             x=0.09,ha="left",fontsize=13,fontweight="bold",y=1.01)
fig.text(0.5,-0.03,"⚠️  There is no correct resolution. Choose it against marker genes and known biology, and report what you chose.",
         ha="center",fontsize=9.5,color="#666")
fig.savefig(OUT+"umap_resolution.png"); print("ok k3")
