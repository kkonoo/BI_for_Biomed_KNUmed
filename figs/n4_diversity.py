from style import *
rng=np.random.default_rng(29)
fig,axes=plt.subplots(1,3,figsize=(13.4,4.4),gridspec_kw={"wspace":0.32})

ax=axes[0]
depth=np.arange(1,40001,200)
for lab,S,c in [("rich sample",320,TEAL),("shallow-sequenced\nsample",320,"#C9A227"),("low-diversity sample",70,RED)]:
    k=0.00016 if "shallow" not in lab else 0.00016
    curve=S*(1-np.exp(-k*depth))
    if "shallow" in lab:
        mask=depth<=8000; ax.plot(depth[mask],curve[mask],color=c,lw=2.4)
        ax.scatter([8000],[curve[mask][-1]],s=55,color=c,zorder=4)
        ax.text(8600,curve[mask][-1],"stopped here",fontsize=9,color=c,va="center")
    else:
        ax.plot(depth,curve,color=c,lw=2.4,label=lab)
ax.set_xlabel("sequencing depth (reads)"); ax.set_ylabel("observed richness")
ax.set_title("① Rarefaction curves",loc="left",fontsize=11.5)
ax.legend(frameon=False,fontsize=9,loc="lower right")
ax.text(0.03,0.96,"Observed richness depends on depth.\nUnequal depth ⇒ unequal richness,\nwith no biology involved.",transform=ax.transAxes,
        va="top",fontsize=9.1,color="#333",linespacing=1.5)

ax=axes[1]
comm=[("A",[0.25]*4),("B",[0.70,0.15,0.10,0.05]),("C",[0.55,0.35,0.06,0.04])]
cols4=[BLUE,TEAL,"#C9A227",PURPLE]
X=np.arange(3); bl=np.zeros(3)
V=np.array([c[1] for c in comm])
for i,c in enumerate(cols4):
    ax.bar(X,V[:,i],bottom=bl,color=c,width=0.6,edgecolor="white",lw=1.1); bl+=V[:,i]
ax.set_xticks(X); ax.set_xticklabels([c[0] for c in comm]); ax.set_ylim(0,1.30)
ax.set_ylabel("relative abundance")
ax.set_title("② Alpha diversity",loc="left",fontsize=11.5)
ax.text(0.0,1.05,"same richness, different evenness",transform=ax.transAxes,fontsize=9.3,color="#666")
for j,v in enumerate(V):
    sh=-(v*np.log(v)).sum()
    ax.text(j,1.05,f"S = 4\nH' = {sh:.2f}",ha="center",fontsize=9.2,color="#444",linespacing=1.4)

ax=axes[2]
g1=rng.multivariate_normal([-1.1,0.3],[[0.55,0.1],[0.1,0.45]],55)
g2=rng.multivariate_normal([1.0,-0.3],[[0.75,-0.15],[-0.15,0.60]],55)
ax.scatter(*g1.T,s=26,color=GREY,alpha=0.85,lw=0,label="Controls")
ax.scatter(*g2.T,s=26,color=RED,alpha=0.7,lw=0,label="Cases")
for g,c in [(g1,"#888"),(g2,RED)]:
    from matplotlib.patches import Ellipse
    m=g.mean(0); cov=np.cov(g.T); w,v=np.linalg.eigh(cov)
    ang=np.degrees(np.arctan2(*v[:,1][::-1]))
    ax.add_patch(Ellipse(m,2*2*np.sqrt(w[1]),2*2*np.sqrt(w[0]),angle=ang,
                 fc="none",ec=c,lw=1.6,ls="--"))
ax.set_xlabel("PCoA 1 (14.2%)"); ax.set_ylabel("PCoA 2 (8.1%)")
ax.set_title("③ Beta diversity",loc="left",fontsize=11.5)
ax.text(0.0,1.05,"PCoA on Bray–Curtis distances",transform=ax.transAxes,fontsize=9.3,color="#666")
ax.legend(frameon=False,fontsize=9.2,loc="upper left")
ax.text(0.985,0.05,"PERMANOVA  R² = 0.06,  p = 0.001\n\nsignificant, and explains 6% of variation —\nreport both numbers",
        transform=ax.transAxes,ha="right",fontsize=9,color="#333",linespacing=1.45,
        bbox=dict(boxstyle="round,pad=0.4",fc="#FAFAFA",ec="#DDD"))
ax.set_xticks([]); ax.set_yticks([])
fig.savefig(OUT+"diversity.png"); print("ok n4")
