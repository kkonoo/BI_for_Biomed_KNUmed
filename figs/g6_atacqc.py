from style import *
rng=np.random.default_rng(23)
fig,axes=plt.subplots(1,2,figsize=(11.5,4.6),gridspec_kw={"wspace":0.26})

x=np.linspace(0,900,900)
def lognorm(c,s,h): return h*np.exp(-0.5*((np.log(np.maximum(x,1))-np.log(c))/s)**2)
sig=lognorm(60,0.45,1.0)+lognorm(215,0.20,0.42)+lognorm(400,0.16,0.16)+lognorm(590,0.13,0.055)
ax=axes[0]
ax.fill_between(x,0,sig,color=TEAL,alpha=0.6,lw=0); ax.plot(x,sig,color=TEAL,lw=1.4)
ax.set_xlabel("fragment length (bp)"); ax.set_ylabel("fragment density")
ax.set_yticks([]); ax.set_xlim(0,900)
ax.set_title("ATAC-seq fragment size distribution",loc="left",fontsize=12)
for c,lab in [(58,"nucleosome-free\n(<100 bp)"),(215,"mono-\nnucleosome"),(400,"di-"),(590,"tri-")]:
    ax.annotate(lab,xy=(c,np.interp(c,x,sig)),xytext=(0,16),textcoords="offset points",
                ha="center",fontsize=9,color="#444",linespacing=1.3)
ax.text(0.985,0.60,"Clear ~200 bp periodicity means the\ntransposition worked.  A smooth,\nfeatureless curve means it did not.",
        transform=ax.transAxes,ha="right",fontsize=9.3,color="#333",linespacing=1.5,
        bbox=dict(boxstyle="round,pad=0.45",fc="#F3F8F8",ec="#BFD8D8"))

ax=axes[1]
d=np.linspace(-2000,2000,800)
good=1+9.5*np.exp(-0.5*(d/230)**2)+rng.normal(0,0.12,800)
bad =1+1.6*np.exp(-0.5*(d/420)**2)+rng.normal(0,0.10,800)
ax.plot(d,good,color=TEAL,lw=2.0,label="good library  (TSS enrichment ≈ 10)")
ax.plot(d,bad,color=RED,lw=2.0,label="failed library  (TSS enrichment ≈ 2.6)")
ax.axvline(0,color="#BBB",lw=1,ls="--")
ax.set_xlabel("distance from transcription start site (bp)")
ax.set_ylabel("normalised signal")
ax.set_title("TSS enrichment — the single best QC number",loc="left",fontsize=12)
ax.legend(frameon=False,fontsize=9.3,loc="upper left")
ax.text(0.985,0.28,"Promoters are open in essentially\nevery cell type, so a flat TSS profile\nmeans the assay failed — not that the\nbiology is unusual.",
        transform=ax.transAxes,ha="right",fontsize=9.3,color="#333",linespacing=1.5,
        bbox=dict(boxstyle="round,pad=0.45",fc="#FAFAFA",ec="#DDD"))
fig.savefig(OUT+"atac_qc.png"); print("ok g6")
