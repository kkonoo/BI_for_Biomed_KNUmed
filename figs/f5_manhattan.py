from style import *
rng=np.random.default_rng(3)
lens=np.array([248,242,198,190,181,171,159,145,138,134,135,133,114,107,102,90,83,80,59,64,47,51],float)
lens/=lens.sum()
nsnp=(lens*400000).astype(int)

fig,ax=plt.subplots(figsize=(12,4.6))
pos0=0; ticks=[]; cols=["#3C6E9C","#9FB6CA"]
peaks={1:[(0.42,11.5)],2:[(0.70,8.9)],5:[(0.30,14.8)],6:[(0.32,31.0)],
       9:[(0.55,9.6)],11:[(0.66,8.3)],16:[(0.24,12.7)],19:[(0.40,10.1)]}
for i,n in enumerate(nsnp):
    chrom=i+1
    x=pos0+np.sort(rng.uniform(0,n,n))
    p=-np.log10(rng.uniform(0,1,n))
    for frac,top in peaks.get(chrom,[]):
        c=int(n*frac); w=max(int(n*0.006),40)
        idx=np.arange(max(0,c-w),min(n,c+w))
        d=np.abs(idx-c)/w
        p[idx]=np.maximum(p[idx],top*np.exp(-2.6*d**1.5)*rng.uniform(0.55,1.0,len(idx)))
    ax.scatter(x,p,s=3.2,color=cols[i%2],lw=0,rasterized=True)
    ticks.append(pos0+n/2); pos0+=n+int(4e5*0.004)

ax.axhline(-np.log10(5e-8),color=RED,lw=1.4,ls="--")
ax.text(pos0*0.999,-np.log10(5e-8)+0.55,"genome-wide significance,  p = 5×10⁻⁸",
        ha="right",fontsize=9.5,color=RED)
ax.axhline(-np.log10(1e-5),color="#C9A227",lw=1.0,ls=":")
ax.text(pos0*0.999,-np.log10(1e-5)+0.35,"suggestive, p = 1×10⁻⁵",ha="right",fontsize=8.5,color="#B08D1F")
ax.set_xticks(ticks); ax.set_xticklabels([str(i+1) for i in range(22)],fontsize=8.5)
ax.set_xlim(-2000,pos0); ax.set_ylim(0,34)
ax.set_xlabel("Chromosome"); ax.set_ylabel(r"$-\log_{10}(p)$")
ax.set_title("Manhattan plot — one point per variant, one column per chromosome",loc="left")
ax.annotate("a locus, not a variant:\nneighbouring SNPs are\ncorrelated by LD",
            xy=(ticks[5],27),xytext=(ticks[8]+9000,29.5),fontsize=9.2,color="#444",
            arrowprops=dict(arrowstyle="-|>",color="#888",lw=1.2,connectionstyle="arc3,rad=-0.2"))
fig.savefig(OUT+"manhattan.png"); print("ok f5")
