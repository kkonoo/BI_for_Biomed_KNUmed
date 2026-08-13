from style import *
rng=np.random.default_rng(12)
fig,axes=plt.subplots(1,2,figsize=(11.8,4.8),gridspec_kw={"wspace":0.28})

# knee plot
nreal=5000; namb=60000
real=rng.lognormal(8.4,0.45,nreal)
amb=rng.lognormal(3.6,0.75,namb)
allc=np.sort(np.concatenate([real,amb]))[::-1]
ax=axes[0]
ax.plot(np.arange(1,len(allc)+1),allc,color=NAVY,lw=2.0)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("barcode rank"); ax.set_ylabel("UMI count per barcode")
ax.set_title("Knee plot — where do the cells stop?",loc="left",fontsize=12)
ax.axvline(nreal,color=RED,ls="--",lw=1.5)
ax.text(nreal*1.25,allc.max()*0.55,"the 'knee'",color=RED,fontsize=10,fontweight="bold")
ax.axvspan(1,nreal,color=TEAL,alpha=0.10,lw=0)
ax.axvspan(nreal,len(allc),color=GREY,alpha=0.18,lw=0)
ax.text(70,allc.max()*0.12,"real cells",color="#2E5E4E",fontsize=10.5,fontweight="bold")
ax.text(2.2e4,allc.max()*0.0025,"empty droplets\ncontaining ambient RNA",color="#666",
        fontsize=9.5,ha="center",linespacing=1.4)
ax.text(0.985,0.055,"⚠️  A hard cutoff here throws away real\nlow-RNA cells. EmptyDrops tests each\nbarcode against the ambient profile instead.",
        transform=ax.transAxes,ha="right",fontsize=9,color="#333",linespacing=1.5,
        bbox=dict(boxstyle="round,pad=0.4",fc="#FAFAFA",ec="#DDD"))

# QC scatter
ax=axes[1]
n=3000
umi=rng.lognormal(8.4,0.42,n)
genes=np.clip(umi**0.62*rng.lognormal(0,0.10,n),0,None)
mito=np.clip(rng.beta(1.6,26,n)*100,0,None)
dying=rng.choice(n,220,replace=False)
mito[dying]=rng.uniform(18,52,220); genes[dying]*=0.55
dbl=rng.choice(np.setdiff1d(np.arange(n),dying),150,replace=False)
umi[dbl]*=2.05; genes[dbl]*=1.55
sc=ax.scatter(umi,genes,c=mito,s=9,cmap="YlOrRd",vmin=0,vmax=40,lw=0,alpha=0.75)
cb=fig.colorbar(sc,ax=ax,pad=0.02); cb.set_label("mitochondrial %",fontsize=9.5)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("UMIs per cell"); ax.set_ylabel("genes detected per cell")
ax.set_title("The standard QC plot",loc="left",fontsize=12)
ax.set_ylim(30,900)
ax.annotate("high mito % =\ndying / lysed cells",
            xy=(float(np.median(umi[dying])),float(np.median(genes[dying]))),
            xytext=(1.35e3,42),fontsize=9.3,color="#B04A42",ha="center",linespacing=1.4,
            arrowprops=dict(arrowstyle="-|>",color="#C1544B",lw=1.2,
                            connectionstyle="arc3,rad=-0.25"))
ax.annotate("high UMI + high genes\n= likely doublets",
            xy=(float(np.median(umi[dbl])),float(np.median(genes[dbl]))),
            xytext=(1.9e3,660),fontsize=9.3,color=NAVY,ha="center",linespacing=1.4,
            arrowprops=dict(arrowstyle="-|>",color=NAVY,lw=1.2,
                            connectionstyle="arc3,rad=0.2"))
fig.savefig(OUT+"knee_qc.png"); print("ok k2")
