from style import *
from scipy.stats import chi2, norm
rng=np.random.default_rng(5)
N=200000
LIM=9.2

def qq(ax,p,title,color,note):
    p=np.sort(p); obs=-np.log10(p); exp=-np.log10((np.arange(1,len(p)+1)-0.5)/len(p))
    lam=np.median(chi2.isf(p,1))/0.4549364
    k=np.unique(np.concatenate([np.arange(0,2500),np.linspace(2500,len(p)-1,4000).astype(int)]))
    ax.plot([0,LIM],[0,LIM],color="#AAAAAA",lw=1.3,ls="--",zorder=1)
    ax.scatter(exp[k],obs[k],s=5,color=color,lw=0,zorder=2)
    ax.set_xlabel(r"Expected $-\log_{10}(p)$"); ax.set_ylabel(r"Observed $-\log_{10}(p)$")
    ax.set_title(title,loc="left",fontsize=12,color=color)
    ax.text(0.05,0.95,f"λ = {lam:.2f}",transform=ax.transAxes,fontsize=14,fontweight="bold",
            color=color,va="top")
    ax.text(0.05,0.855,note,transform=ax.transAxes,fontsize=9.4,color="#444",va="top",linespacing=1.55)
    ax.set_xlim(0,6.0); ax.set_ylim(0,LIM)

fig,axes=plt.subplots(1,2,figsize=(11.5,5.2),gridspec_kw={"wspace":0.26})

p_good=rng.uniform(0,1,N)
hit=rng.choice(N,140,replace=False)
p_good[hit]=10**(-rng.gamma(3.0,1.6,140))
qq(axes[0],p_good,"Well-controlled study",BLUE,
   "points sit on the diagonal until the\nvery tail — that small departure\nat the top is the real signal")

z=rng.normal(0,1.19,N); p_bad=2*norm.sf(np.abs(z))
qq(axes[1],p_bad,"⚠️  Inflated — something is wrong",RED,
   "the entire distribution lifts off the\ndiagonal: unmodelled population\nstructure, cryptic relatedness,\nor batch effects")

fig.suptitle("QQ plots: the first figure to look at in any GWAS",x=0.09,ha="left",
             fontsize=13,fontweight="bold",y=1.02)
fig.text(0.5,-0.02,"λ (genomic inflation factor) ≈ 1 means calibrated. λ > 1.1 means the p-values cannot be trusted as they stand.",
         ha="center",fontsize=9,color="#666")
fig.savefig(OUT+"qq_plot.png"); print("ok f6")
