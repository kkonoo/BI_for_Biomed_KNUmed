from style import *
rng=np.random.default_rng(11)
fig,axes=plt.subplots(1,2,figsize=(11.5,5.0),gridspec_kw={"wspace":0.28})

pops=[("AFR",(-2.6,0.15),(0.62,0.55),"#C1544B"),
      ("EUR",(1.55,1.35),(0.30,0.28),"#4A77A8"),
      ("EAS",(1.75,-1.25),(0.26,0.26),"#5F9EA0"),
      ("SAS",(1.05,0.05),(0.34,0.30),"#B5892A"),
      ("AMR",(0.15,0.85),(0.55,0.45),"#9B8EC4")]
ax=axes[0]
for lab,(mx,my),(sx,sy),c in pops:
    n=260
    ax.scatter(rng.normal(mx,sx,n),rng.normal(my,sy,n),s=9,color=c,alpha=0.55,lw=0,label=lab)
ax.set_xlabel("PC1  (18.2% of variance)"); ax.set_ylabel("PC2  (6.4%)")
ax.set_title("Genotype PCA recovers continental ancestry",loc="left",fontsize=12)
ax.legend(frameon=False,fontsize=9.5,markerscale=1.9,loc="lower left",ncol=2,
          handletextpad=0.2,columnspacing=0.9)
ax.set_xticks([]); ax.set_yticks([])

ax=axes[1]
n=300
eur=rng.normal(1.55,0.30,n); afr=rng.normal(-2.6,0.62,n)
cases=np.concatenate([rng.normal(1.55,0.30,90),rng.normal(-2.6,0.62,210)])
ctrl =np.concatenate([rng.normal(1.55,0.30,220),rng.normal(-2.6,0.62,80)])
ax.hist(ctrl,bins=34,color=GREY,alpha=0.75,label="Controls",lw=0)
ax.hist(cases,bins=34,color=RED,alpha=0.62,label="Cases",lw=0)
ax.set_xlabel("PC1"); ax.set_ylabel("number of individuals")
ax.set_title("⚠️  When ancestry differs between groups",loc="left",fontsize=12,color=RED)
ax.legend(frameon=False,fontsize=9.5,loc="upper center")
ax.text(0.5,0.60,"Any variant whose frequency differs\nbetween these ancestries will appear\nassociated with disease — with no\ncausal relationship whatsoever.",
        transform=ax.transAxes,ha="center",va="top",fontsize=9.6,color="#444",linespacing=1.5,
        bbox=dict(boxstyle="round,pad=0.5",fc="#FDF3F2",ec="#E8C4C0"))
ax.set_yticks([])
fig.savefig(OUT+"pca_population.png"); print("ok f4")
