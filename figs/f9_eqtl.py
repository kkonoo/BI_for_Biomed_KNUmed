from gallery import *
rng=np.random.default_rng(9)
fig,ax=plt.subplots(figsize=(8.6,4.6))
means=[3.2,5.6,8.1]
for g,(lab,m) in enumerate(zip(["AA","AB","BB"],means)):
    v=rng.normal(m,1.1,60)
    ax.scatter(g+rng.uniform(-0.18,0.18,60),v,s=13,color="#C9CED6",lw=0)
    ax.plot([g-0.25,g+0.25],[m,m],color=CRIMSON,lw=2.6)
ax.plot([0,1,2],means,color=CRIMSON,lw=1,ls=(0,(3,2)))
ax.text(2.33,means[2],"mean",va="center",fontsize=9.5,color=CRIMSON)
ax.text(2.33,5.0,"one dot = one person\n\nthe more B alleles,\nthe more gene G is made\n→ the variant is an eQTL for G",
        va="center",fontsize=9.8,color=TEXT2,linespacing=1.35)
ax.set_xlim(-0.5,3.4); ax.set_xticks([0,1,2],["AA","AB","BB"]); ax.set_yticks([])
ax.spines["bottom"].set_bounds(-0.5,2.5)
for s in ("left","bottom"): ax.spines[s].set_color(LINE)
ax.tick_params(colors=TEXT2)
ax.set_xlabel("genotype at the GWAS variant",color=TEXT2); ax.set_ylabel("expression of gene G",color=TEXT2)
header(ax,"eQTL","Does the variant change how much of a gene is made?")
fig.savefig(OUT+"eqtl.png"); print("ok eqtl")
