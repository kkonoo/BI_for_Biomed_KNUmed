from style import *
import matplotlib.patches as mp
rng=np.random.default_rng(41)
N=220
tads=[(0,58),(58,104),(104,168),(168,220)]
M=np.zeros((N,N))
ii,jj=np.meshgrid(np.arange(N),np.arange(N),indexing="ij")
M=1.0/(1.0+np.abs(ii-jj)/3.0)**1.05
for a,b in tads:
    blk=(ii>=a)&(ii<b)&(jj>=a)&(jj<b)
    M[blk]*=2.5
M+= 0.02*rng.random((N,N)); M=(M+M.T)/2
# a strong loop between an enhancer and a promoter
for (a,b,s) in [(118,160,0.65)]:
    for da in range(-3,4):
        for db in range(-3,4):
            w=np.exp(-(da*da+db*db)/6)
            M[a+da,b+db]+=s*w; M[b+db,a+da]+=s*w

fig,(ax,gx)=plt.subplots(2,1,figsize=(7.8,7.6),
                         gridspec_kw={"height_ratios":[5,0.85],"hspace":0.06})
im=ax.imshow(np.log1p(M*6),cmap="Reds",origin="upper",interpolation="nearest",aspect="auto")
for a,b in tads:
    ax.plot([a,b,a,a],[a,b,b,a],color=NAVY,lw=1.5,ls="--")
ax.plot([160],[118],marker="o",ms=15,mfc="none",mec="#1F5C8B",mew=2.2)
ax.annotate("enhancer–promoter loop",xy=(165,112),xytext=(112,52),fontsize=9.5,color="#1F5C8B",
            arrowprops=dict(arrowstyle="-|>",color="#1F5C8B",lw=1.3))
ax.text(42,16,"TAD",fontsize=10,color=NAVY,ha="center",va="center",fontweight="bold")
ax.text(0.98,0.06,"contacts are frequent within a TAD\nand rare across its boundary",
        transform=ax.transAxes,ha="right",fontsize=9.3,color="#333",linespacing=1.5,
        bbox=dict(boxstyle="round,pad=0.45",fc="white",ec="#DDD",alpha=0.92))
ax.set_xticks([]); ax.set_yticks([])
ax.set_title("Hi-C contact matrix",loc="left",fontsize=12.5)
for sp in ax.spines.values(): sp.set_visible(True); sp.set_color("#AAA")

gx.plot([0,N],[0.62,0.62],color="#9AA5B1",lw=1.6)
for s,w,c,lab in [(20,26,NAVY,"GENE-B"),(112,10,"#C9A227","enhancer"),(155,22,"#1F5C8B","GENE-C")]:
    gx.add_patch(mp.Rectangle((s,0.50),w,0.24,fc=c,ec="none"))
    gx.text(s+w/2,0.30,lab,ha="center",fontsize=8.8,color=c)
for a,b in tads: gx.axvline(b,color=NAVY,lw=1.2,ls="--",alpha=0.6)
gx.text(N/2,-0.16,"TAD boundaries (dashed) are largely shared between cell types; loops inside them are not",
        ha="center",fontsize=8.8,color="#666")
gx.set_xlim(0,N); gx.set_ylim(0,1); gx.axis("off")
fig.savefig(OUT+"hic_tad.png"); print("ok g7")
