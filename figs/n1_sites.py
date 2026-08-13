from style import *
rng=np.random.default_rng(9)
fig,axes=plt.subplots(1,2,figsize=(12.6,4.8),gridspec_kw={"wspace":0.30,"width_ratios":[1.15,1]})
taxa=["Bacteroidetes","Firmicutes","Proteobacteria","Actinobacteria","Fusobacteria","Other"]
cols=[BLUE,TEAL,"#C9A227","#9B8EC4","#E08A3C","#CFCFCF"]
sites={"Gut":[0.44,0.42,0.05,0.04,0.01,0.04],
       "Oral":[0.16,0.24,0.14,0.12,0.22,0.12],
       "Skin":[0.06,0.24,0.16,0.46,0.01,0.07],
       "Vaginal":[0.02,0.90,0.02,0.02,0.01,0.03],
       "Airway":[0.14,0.36,0.28,0.10,0.05,0.07]}
ax=axes[0]
X=np.arange(len(sites)); bl=np.zeros(len(sites))
V=np.array(list(sites.values()))
for i,(t,c) in enumerate(zip(taxa,cols)):
    ax.bar(X,V[:,i],bottom=bl,color=c,width=0.62,label=t,edgecolor="white",lw=1.1); bl+=V[:,i]
ax.set_xticks(X); ax.set_xticklabels(sites.keys())
ax.set_ylabel("relative abundance (phylum)"); ax.set_ylim(0,1)
ax.set_title("Body sites host entirely different communities",loc="left",fontsize=12)
ax.legend(frameon=False,fontsize=8.8,loc="upper center",bbox_to_anchor=(0.5,-0.13),ncol=3,
          handletextpad=0.35,columnspacing=0.9)

ax=axes[1]
n=14
base=np.array([0.40,0.44,0.06,0.05,0.01,0.04])
M=np.array([rng.dirichlet(base*7.0) for _ in range(n)])
bl=np.zeros(n)
for i,c in enumerate(cols):
    ax.bar(np.arange(n),M[:,i],bottom=bl,color=c,width=0.82,edgecolor="white",lw=0.7); bl+=M[:,i]
ax.set_xticks([]); ax.set_ylim(0,1); ax.set_ylabel("relative abundance")
ax.set_xlabel("14 healthy adults, gut")
ax.set_title("…and healthy people differ enormously from each other",loc="left",fontsize=12)
ax.text(0.5,-0.28,"⚠️  With variation this large between healthy individuals, a small case–control study has very little power.",
        transform=ax.transAxes,ha="center",va="top",fontsize=9.5,color="#8A3A33",linespacing=1.4)
fig.savefig(OUT+"taxonomic_composition.png"); print("ok n1")
