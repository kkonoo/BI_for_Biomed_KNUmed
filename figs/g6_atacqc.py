from gallery import *
rng=np.random.default_rng(23)
fig,axes=plt.subplots(1,2,figsize=(11.5,4.6),gridspec_kw={"wspace":0.22})

# (a) fragment sizes: nucleosome-free peak, then ~200 bp steps (mono-, di-, tri-nucleosome)
x=np.linspace(0,900,900)
def lognorm(c,s,h): return h*np.exp(-0.5*((np.log(np.maximum(x,1))-np.log(c))/s)**2)
sig=lognorm(60,0.45,1.0)+lognorm(215,0.20,0.42)+lognorm(400,0.16,0.16)+lognorm(590,0.13,0.055)
ax=axes[0]
ax.fill_between(x,0,sig,color=tint(0.25),lw=0); ax.plot(x,sig,color=CRIMSON,lw=1.6)
ax.set_xlabel("fragment length (bp)",color=TEXT2); ax.set_ylabel("fragment density",color=TEXT2)
ax.set_yticks([]); ax.set_xlim(0,900); ax.set_ylim(0,1.22)
ax.text(88,1.0,"nucleosome-free  (< 100 bp)",ha="left",va="center",fontsize=9.5,color=TEXT)
for c,lab,dy in [(215,"mono-\nnucleosome",10),(400,"di-",8),(590,"tri-",8)]:
    ax.annotate(lab,xy=(c,np.interp(c,x,sig)),xytext=(0,dy),textcoords="offset points",
                ha="center",va="bottom",fontsize=9.5,color=TEXT,linespacing=1.25)
ax.text(0.98,0.56,"~200 bp steps → transposition worked\nsmooth, featureless curve → it failed",
        transform=ax.transAxes,ha="right",va="center",fontsize=9.5,color=TEXT2,linespacing=1.5)
ax.set_title("Fragment size distribution",loc="left",fontsize=12.5,fontweight="semibold",color=TEXT,pad=8)

# (b) TSS enrichment: good ≈ 10, failed ≈ 2.6
bx=axes[1]
d=np.linspace(-2000,2000,800)
good=1+9.5*np.exp(-0.5*(d/230)**2)+rng.normal(0,0.12,800)
bad =1+1.6*np.exp(-0.5*(d/420)**2)+rng.normal(0,0.10,800)
bx.axvline(0,color=LINE,lw=1,ls=(0,(3,3)))
bx.plot(d,bad,color=TEXT2,lw=1.6)
bx.plot(d,good,color=CRIMSON,lw=1.8)
bx.text(260,9.6,"good library\nTSS enrichment ≈ 10",ha="left",va="center",fontsize=9.5,color=CRIMSON,
        fontweight="semibold",linespacing=1.3)
bx.text(560,3.15,"failed library\nTSS enrichment ≈ 2.6",ha="left",va="bottom",fontsize=9.5,color=TEXT2,linespacing=1.3)
bx.text(0.02,0.62,"promoters are open in\nalmost every cell type\n→ flat = assay failure",
        transform=bx.transAxes,ha="left",va="top",fontsize=9.5,color=TEXT2,linespacing=1.5)
bx.set_xticks([-2000,-1000,0,1000,2000])
bx.set_xlabel("distance from the transcription start site (bp)",color=TEXT2)
bx.set_ylabel("normalised signal",color=TEXT2)
bx.set_xlim(-2000,2000); bx.set_ylim(0,11.5)
bx.set_title("TSS enrichment",loc="left",fontsize=12.5,fontweight="semibold",color=TEXT,pad=8)

for a in axes:
    for s in ("left","bottom"): a.spines[s].set_color(LINE)
    a.tick_params(colors=TEXT2)
header(axes[0],"ATAC-seq QC","A working library shows nucleosome steps and a sharp peak at promoters",dy=24,x1=2.2)
fig.savefig(OUT+"atac_qc.png"); print("ok g6")
