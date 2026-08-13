from style import *
rng=np.random.default_rng(55)
def blob(n,c,s,rot=0):
    p=np.stack([rng.normal(0,s[0],n),rng.normal(0,s[1],n)],1)
    co,si=np.cos(rot),np.sin(rot); p=p@np.array([[co,-si],[si,co]])
    return p+np.array(c)
cts=[("T",(-2.6,2.4),(1.2,0.8),0.4),("B",(2.8,1.4),(0.85,0.7),0),
     ("Mono",(1.6,-2.8),(1.1,0.9),-0.35),("NK",(-1.0,-2.0),(0.6,0.5),0)]
ctcol={"T":BLUE,"B":TEAL,"Mono":"#C9A227","NK":PURPLE}
shift={"A":np.array([0,0]),"B":np.array([2.6,2.2])}
def build(sep):
    P=[];CT=[];BA=[]
    for b,off in shift.items():
        for nm,c,s,r in cts:
            n=500
            P.append(blob(n,np.array(c)+off*sep,s,r)); CT+= [nm]*n; BA+=[b]*n
    return np.vstack(P),np.array(CT),np.array(BA)

fig,axes=plt.subplots(2,3,figsize=(13.4,7.6),gridspec_kw={"wspace":0.10,"hspace":0.28})
states=[("① Before integration",1.0),("② Integrated",0.0),("③ ⚠️ Over-corrected",0.0)]
for j,(t,sep) in enumerate(states):
    P,CT,BA=build(sep)
    if j==2:
        P=P*np.array([0.55,0.55])+rng.normal(0,0.85,P.shape)
    for i,key in enumerate(["batch","celltype"]):
        ax=axes[i,j]
        if key=="batch":
            for b,c in [("A","#4A77A8"),("B","#E08A3C")]:
                m=BA==b; ax.scatter(P[m,0],P[m,1],s=4,color=c,lw=0,alpha=0.6,label=f"donor {b}")
            if j==0: ax.legend(frameon=False,fontsize=9,markerscale=2.4,loc="lower left")
        else:
            for nm in ctcol:
                m=CT==nm; ax.scatter(P[m,0],P[m,1],s=4,color=ctcol[nm],lw=0,alpha=0.65,label=nm)
            if j==0: ax.legend(frameon=False,fontsize=9,markerscale=2.4,loc="lower left",ncol=2,
                               handletextpad=0.2,columnspacing=0.7)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlim(-8,9); ax.set_ylim(-7,8)
        if i==0: ax.set_title(t,loc="left",fontsize=12,color=RED if j==2 else "#222")
        if j==0: ax.set_ylabel("coloured by donor" if key=="batch" else "coloured by cell type",fontsize=10.5)
notes=["cells separate by donor before\nthey separate by cell type",
       "donors mixed, cell types preserved\n— this is what you want",
       "donors mixed, but the cell types\nare gone too"]
for j,(nt,c) in enumerate(zip(notes,["#555","#2E5E4E","#8A3A33"])):
    axes[1,j].text(0.5,-0.13,nt,transform=axes[1,j].transAxes,ha="center",va="top",
                   fontsize=9.5,color=c,linespacing=1.4)
fig.text(0.5,-0.02,"⚠️  Always inspect both colourings. A plot coloured only by donor cannot tell you whether integration destroyed the biology.",
         ha="center",fontsize=9.8,color="#666")
fig.savefig(OUT+"batch_integration.png"); print("ok k4")
