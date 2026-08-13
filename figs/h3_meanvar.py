from style import *
rng=np.random.default_rng(8)
fig,axes=plt.subplots(1,2,figsize=(11.5,4.7),gridspec_kw={"wspace":0.26})
G=4000
mu=np.exp(rng.normal(3.4,1.9,G))
disp=0.02+18/(mu+30)
var=mu+disp*mu**2
obs=var*np.exp(rng.normal(0,0.28,G))
ax=axes[0]
ax.scatter(mu,obs,s=5,color=BLUE,alpha=0.30,lw=0)
m=np.logspace(0,5,200)
ax.plot(m,m,color="#666",lw=2.0,ls="--",label="Poisson:  variance = mean")
ax.plot(m,m+(0.02+18/(m+30))*m**2,color=RED,lw=2.2,label="Negative binomial:  var = μ + φμ²")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("mean count across replicates"); ax.set_ylabel("variance across replicates")
ax.set_title("Counts are overdispersed",loc="left",fontsize=12)
ax.legend(frameon=False,fontsize=9.2,loc="upper left")
ax.text(0.985,0.06,"Real replicates vary far more than\nPoisson allows. Using a Poisson or a\nt-test on counts under-estimates\nvariance → false positives.",
        transform=ax.transAxes,ha="right",fontsize=9.2,color="#333",linespacing=1.5,
        bbox=dict(boxstyle="round,pad=0.45",fc="#FAFAFA",ec="#DDD"))

ax=axes[1]
raw=disp*np.exp(rng.normal(0,0.55,G))
trend=0.02+18/(mu+30)
shrunk=np.exp(0.35*np.log(raw)+0.65*np.log(trend))
o=np.argsort(mu)
ax.scatter(mu,raw,s=6,color="#C9C9C9",alpha=0.55,lw=0,label="gene-wise estimate (noisy)")
ax.plot(mu[o],trend[o],color=NAVY,lw=2.4,label="fitted trend")
ax.scatter(mu,shrunk,s=6,color=TEAL,alpha=0.55,lw=0,label="shrunken estimate (used)")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("mean count"); ax.set_ylabel("dispersion  φ")
ax.set_title("Why DESeq2 / edgeR borrow information",loc="left",fontsize=12)
ax.legend(frameon=False,fontsize=9.2,loc="lower left")
ax.text(0.985,0.94,"With n = 3 you cannot estimate a\ngene's variance. These methods\nshrink each gene toward the trend\nfitted from all genes.",
        transform=ax.transAxes,ha="right",va="top",fontsize=9.2,color="#333",linespacing=1.5,
        bbox=dict(boxstyle="round,pad=0.45",fc="#F3F8F8",ec="#BFD8D8"))
fig.savefig(OUT+"mean_variance.png"); print("ok h3")
