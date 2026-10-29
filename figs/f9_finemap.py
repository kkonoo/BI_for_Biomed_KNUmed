from gallery import *
# same simulated locus as f7_locus.py (regional association plot)
rng=np.random.default_rng(21)
lo,hi=44.0,44.9
lead=44.46
n=1400
pos=np.sort(rng.uniform(lo,hi,n))
d=np.abs(pos-lead)
r2=np.exp(-((d/0.055)**1.7))*rng.uniform(0.75,1.0,n)
r2=np.clip(r2,0,1)
logp=r2*13.5+rng.gamma(1.0,0.85,n)
logp=np.clip(logp,0.02,None)
i=np.argmin(np.abs(pos-lead)); logp[i]=logp.max()+0.7; r2[i]=1.0

# fine-mapping with one causal variant: posterior ∝ approximate Bayes factor (Wakefield)
from scipy.stats import norm
z=norm.isf(10**-logp/2)
V,W=0.02**2,0.15**2                      # sampling variance of β, prior variance of β
lbf=0.5*np.log(V/(V+W))+z**2/2*W/(V+W)
pip=np.exp(lbf-lbf.max()); pip/=pip.sum()
order=np.argsort(-pip); cs=order[:np.searchsorted(np.cumsum(pip[order]),0.95)+1]
inset=np.zeros(n,bool); inset[cs]=True

fig,(ax,bx)=plt.subplots(2,1,figsize=(10.5,6.6),gridspec_kw={"height_ratios":[1.5,1],"hspace":0.42})
zl,zh=44.43,44.49                       # zoom window for the lower panel
ax.scatter(pos[~inset],logp[~inset],s=18,color="#C9CED6",lw=0)
ax.scatter(pos[inset],logp[inset],s=46,color=CRIMSON,lw=0.6,edgecolor="white",zorder=3)
ax.scatter([pos[i]],[logp[i]],s=120,marker="D",facecolor="white",edgecolor=TEXT,lw=1.3,zorder=4)
ax.annotate("lead SNP",xy=(pos[i],logp[i]),xytext=(pos[i]-0.13,logp[i]+0.3),fontsize=10,color=TEXT,
            fontweight="semibold",arrowprops=dict(arrowstyle="-",color=TEXT2,lw=0.8))
ax.axhline(-np.log10(5e-8),color=TEXT2,lw=0.9,ls=(0,(4,3)))
ax.text(lo+0.005,-np.log10(5e-8)+0.3,"p = 5×10⁻⁸",fontsize=9,color=TEXT2)
ax.set_ylabel(r"$-\log_{10}(p)$",color=TEXT2); ax.set_ylim(0,17.4)
ax.text(0.99,0.95,f"{(logp>-np.log10(5e-8)).sum()} variants pass 5×10⁻⁸",transform=ax.transAxes,
        ha="right",va="top",fontsize=10,color=TEXT2)

ax.axvspan(zl,zh,color=FAINT,zorder=0)
ax.text((zl+zh)/2,0.4,"zoomed below",ha="center",fontsize=8.5,color=TEXT2)
w=(pos>=zl)&(pos<=zh)
bx.vlines(pos[w&~inset],0,np.maximum(pip[w&~inset],0.004),color="#C9CED6",lw=1)
for k in cs:
    bx.plot([pos[k],pos[k]],[0,pip[k]],color=CRIMSON,lw=2)
    bx.scatter([pos[k]],[pip[k]],s=30,color=CRIMSON,zorder=3)
for k in cs:
    if k==i: continue
    left=any(0<pos[m]-pos[k]<0.003 for m in cs)   # a close neighbour on the right → label to the left
    bx.text(pos[k]+(-0.0006 if left else 0.0006),pip[k]+0.035,f"{pip[k]:.2f}",ha="right" if left else "left",
            fontsize=8.5,color=CRIMSON)
bx.text(pos[i]+0.0015,pip[i],f"  {pip[i]:.2f}",ha="left",va="center",fontsize=8.5,color=TEXT)
bx.scatter([pos[i]],[pip[i]],s=90,marker="D",facecolor="white",edgecolor=TEXT,lw=1.2,zorder=4)
bx.set_ylabel("posterior probability\nof being causal",color=TEXT2); bx.set_ylim(0,max(pip)*1.25)
bx.text(0.99,0.92,f"95% credible set = {len(cs)} variants\n(probabilities add up to ≥ 0.95)",transform=bx.transAxes,ha="right",va="top",
        fontsize=11,color=CRIMSON,fontweight="bold")
bx.set_xlabel("position on chromosome 6 (Mb)",color=TEXT2); bx.set_xlim(zl,zh); ax.set_xlim(lo,hi)
for a in (ax,bx):
    for s in ("left","bottom"): a.spines[s].set_color(LINE)
    a.tick_params(colors=TEXT2)
header(ax,"Fine-mapping","From many associated variants to a short list that holds the causal one")
fig.savefig(OUT+"finemap_credible.png"); print("ok f9",len(cs),pip[cs].round(3),pos[cs].round(3),pip[i].round(3))
