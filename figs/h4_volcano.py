from style import *
rng=np.random.default_rng(15)
G=12000
mu=np.exp(rng.normal(3.6,2.0,G))
lfc=rng.normal(0,0.32,G)
idx=rng.choice(G,420,replace=False)
lfc[idx]+=rng.choice([-1,1],420)*rng.gamma(2.4,0.55,420)
se=0.28+2.4/np.sqrt(mu+8)
z=lfc/se
p=2*np.exp(-0.5*z**2)/ (1+0.0)   # rough
p=np.clip(np.exp(-0.5*z**2)*0.9,1e-40,1)
padj=np.clip(p*G/np.maximum(np.argsort(np.argsort(p))+1,1),1e-40,1)
sig=(padj<0.05)&(np.abs(lfc)>1)

fig,axes=plt.subplots(1,2,figsize=(11.8,5.0),gridspec_kw={"wspace":0.26})
ax=axes[0]
col=np.where(sig&(lfc>0),RED,np.where(sig&(lfc<0),BLUE,"#C8C8C8"))
ax.scatter(lfc,-np.log10(padj),s=6,c=col,alpha=0.6,lw=0)
ax.axhline(-np.log10(0.05),color="#888",ls="--",lw=1.1)
ax.axvline(1,color="#888",ls=":",lw=1.0); ax.axvline(-1,color="#888",ls=":",lw=1.0)
ax.set_xlabel("log₂ fold change"); ax.set_ylabel("−log₁₀(adjusted p)")
ax.set_title("Volcano plot",loc="left",fontsize=12)
ax.set_xlim(-6,6)
ax.text(0.02,0.97,f"{(sig&(lfc<0)).sum()} down",transform=ax.transAxes,color=BLUE,fontsize=10,
        fontweight="bold",va="top")
ax.text(0.98,0.97,f"{(sig&(lfc>0)).sum()} up",transform=ax.transAxes,color=RED,fontsize=10,
        fontweight="bold",va="top",ha="right")
ax.text(0.5,0.03,"⚠️  Both dashed lines are conventions, not\nfacts. |log₂FC| > 1 and padj < 0.05 are\nhabits — state and justify yours.",
        transform=ax.transAxes,ha="center",fontsize=9.2,color="#8A3A33",linespacing=1.5,
        bbox=dict(boxstyle="round,pad=0.42",fc="#FDF3F2",ec="#E8C4C0"))

ax=axes[1]
ax.scatter(mu,lfc,s=6,c=col,alpha=0.6,lw=0)
ax.set_xscale("log"); ax.axhline(0,color="#888",lw=1.1)
ax.set_xlabel("mean normalised count"); ax.set_ylabel("log₂ fold change")
ax.set_title("MA plot",loc="left",fontsize=12)
ax.set_ylim(-6,6)
ax.annotate("low-expression genes have huge\napparent fold changes and no power",
            xy=(1.2,-3.6),xytext=(120,-5.2),fontsize=9.2,color="#444",ha="center",linespacing=1.4,
            arrowprops=dict(arrowstyle="-|>",color="#999",lw=1.2,connectionstyle="arc3,rad=0.25"))
ax.text(0.985,0.96,"The MA plot is the better diagnostic:\nit shows whether your 'hits' are just\nnoise from lowly expressed genes.",
        transform=ax.transAxes,ha="right",va="top",fontsize=9.2,color="#333",linespacing=1.5,
        bbox=dict(boxstyle="round,pad=0.42",fc="#FAFAFA",ec="#DDD"))
fig.savefig(OUT+"volcano_ma.png"); print("ok h4")
