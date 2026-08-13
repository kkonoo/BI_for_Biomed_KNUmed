from style import *
rng=np.random.default_rng(19)
fig,axes=plt.subplots(1,2,figsize=(11.8,4.7),gridspec_kw={"wspace":0.28})
ax=axes[0]
biomass=np.logspace(1,9,200)
contam=100/(1+biomass/3e3)
ax.plot(biomass,contam,color=RED,lw=2.6)
ax.set_xscale("log")
ax.set_xlabel("microbial biomass in sample (cells)"); ax.set_ylabel("% of reads from reagent contaminants")
ax.set_title("The kit-ome",loc="left",fontsize=12,color=RED)
ax.set_ylim(-3,103)
zones=[(3e7,"Stool","#2E8B57",6,22),(3e5,"Saliva,\nvaginal swab",TEAL,-4,22),
       (2e3,"Skin swab, BAL","#C9A227",14,10),(60,"Blood, tissue biopsy,\nplacenta, tumour",RED,16,-34)]
for x,lab,c,dx,dy in zones:
    y=100/(1+x/3e3)
    ax.scatter([x],[y],s=70,color=c,zorder=4,ec="white",lw=1.3)
    ax.annotate(lab,(x,y),xytext=(dx,dy),textcoords="offset points",fontsize=9.3,color=c,
                fontweight="bold",linespacing=1.35,ha="left" if dx>=0 else "center")
ax.axhspan(50,103,color=RED,alpha=0.07,lw=0)
ax.text(1.2e6,78,"below ~10³ cells, most of what you\nsequence came out of the kit",
        fontsize=9.4,color="#8A3A33",ha="center",linespacing=1.45)

ax=axes[1]; ax.axis("off")
ax.set_title("What to do about it",loc="left",fontsize=12)
items=[("Run blank controls","Extraction blanks and no-template PCR controls,\nin every batch. Sequence them.",TEAL),
       ("Use them, don't discard them","decontam and similar tools identify contaminants\nby (i) prevalence in blanks and (ii) inverse\ncorrelation with DNA concentration.",TEAL),
       ("Record DNA concentration","It is the key covariate. Contaminant fraction\nscales inversely with input DNA.",TEAL),
       ("⚠️ Distrust low-biomass claims","'The healthy placenta has a microbiome' and\nseveral tumour-microbiome results were later\nattributed largely to contamination.",RED)]
y=0.94
for t,d,c in items:
    ax.text(0.0,y,t,fontsize=10.6,fontweight="bold",color=c,va="top")
    ax.text(0.02,y-0.062,d,fontsize=9.3,color="#555",va="top",linespacing=1.5)
    y-=0.255
ax.set_xlim(0,1); ax.set_ylim(-0.05,1)
fig.savefig(OUT+"contamination.png"); print("ok n2")
