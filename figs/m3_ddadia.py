from style import *
import matplotlib.patches as mp
rng=np.random.default_rng(21)
fig,axes=plt.subplots(2,2,figsize=(12.4,6.6),gridspec_kw={"hspace":0.78,"wspace":0.22,
                                                          "height_ratios":[1,0.85]})
mz=np.linspace(400,1000,900)
peaks=[(438,0.95),(471,0.30),(505,0.62),(552,0.18),(588,1.00),(631,0.24),(668,0.44),
       (702,0.13),(745,0.71),(788,0.21),(824,0.35),(869,0.10),(902,0.52),(946,0.16)]
spec=sum(h*np.exp(-0.5*((mz-c)/1.1)**2) for c,h in peaks)

ax=axes[0,0]
ax.fill_between(mz,0,spec,color="#B9C4CF",lw=0)
top=sorted(peaks,key=lambda t:-t[1])[:4]
for c,h in top:
    ax.vlines(c,0,h,color=RED,lw=2.2)
    ax.plot([c],[h+0.06],marker="v",color=RED,ms=8)
ax.set_ylim(0,1.28); ax.set_xlim(400,1000); ax.set_yticks([])
ax.set_title("DDA — data-dependent acquisition",loc="left",fontsize=12,color=RED)
ax.set_xlabel("MS1  m/z")
ax.text(0.985,0.95,"the instrument picks the top-N most\nintense precursors and fragments those",
        transform=ax.transAxes,ha="right",va="top",fontsize=9.3,color="#333",linespacing=1.5)

ax=axes[0,1]
ax.fill_between(mz,0,spec,color="#B9C4CF",lw=0)
for i,x0 in enumerate(np.arange(400,1000,50)):
    ax.add_patch(mp.Rectangle((x0,0),50,1.28,fc=TEAL,alpha=0.10 if i%2 else 0.20,ec=TEAL,lw=0.8))
ax.set_ylim(0,1.28); ax.set_xlim(400,1000); ax.set_yticks([])
ax.set_title("DIA — data-independent acquisition",loc="left",fontsize=12,color=TEAL)
ax.set_xlabel("MS1  m/z")
ax.text(0.985,0.95,"fixed isolation windows march across the\nwhole range, fragmenting everything",
        transform=ax.transAxes,ha="right",va="top",fontsize=9.3,color="#333",linespacing=1.5)

# consequence: run-to-run missingness
for ax,lab,c,miss in [(axes[1,0],"DDA",RED,0.42),(axes[1,1],"DIA",TEAL,0.11)]:
    R,P=10,40
    M=(rng.random((P,R))<miss)
    ax.imshow(~M,cmap="Greys_r",aspect="auto",vmin=0,vmax=1)
    ax.imshow(np.ma.masked_where(~M,M),cmap=plt.matplotlib.colors.ListedColormap([c]),
              aspect="auto",alpha=0.85)
    ax.set_xlabel("run"); ax.set_ylabel("protein")
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"{lab}:  {miss*100:.0f}% missing across runs",loc="left",fontsize=11.5,color=c)
fig.text(0.5,0.445,"Consequence — how often a protein is measured in every run     (coloured = missing)",
         ha="center",fontsize=11,fontweight="bold",color="#333")
fig.text(0.5,-0.02,"DDA's selection is stochastic, so a low-abundance peptide picked in run 1 may be skipped in run 2. DIA fragments everything, so missingness is far lower — at the cost of much more complex spectra to deconvolute.",
         ha="center",fontsize=9.4,color="#666")
fig.savefig(OUT+"dda_vs_dia.png"); print("ok m3")
