from style import *
rng=np.random.default_rng(19)
N=2000
fig=plt.figure(figsize=(12.6,5.4))
gs=fig.add_gridspec(3,2,height_ratios=[2.6,0.42,1.5],width_ratios=[1.25,1],hspace=0.16,wspace=0.28)

# --- GSEA running score ---
rank=np.arange(N)
inset=np.sort(np.concatenate([rng.choice(np.arange(0,320),34,replace=False),
                              rng.choice(np.arange(320,N),16,replace=False)]))
hit=np.zeros(N,bool); hit[inset]=True
w=np.linspace(2.6,-2.6,N)
Nh=np.abs(w[hit]).sum()
step=np.where(hit,np.abs(w)/Nh,-1/(N-hit.sum()))
run=np.cumsum(step)

ax=fig.add_subplot(gs[0,0])
ax.plot(rank,run,color=TEAL,lw=2.2)
ax.axhline(0,color="#BBB",lw=1)
k=np.argmax(run)
ax.plot([k,k],[0,run[k]],color=RED,lw=1.6,ls="--")
ax.annotate(f"ES = {run[k]:.2f}",xy=(k,run[k]),xytext=(18,-14),textcoords="offset points",
            fontsize=10,fontweight="bold",color=RED)
ax.set_ylabel("running enrichment score")
ax.set_xlim(0,N); ax.set_xticks([])
ax.set_title("GSEA — uses the whole ranked list",loc="left",fontsize=12)
ax.text(0.985,0.06,"no threshold anywhere:\nevery gene contributes",transform=ax.transAxes,
        ha="right",fontsize=9.2,color="#555",style="italic",linespacing=1.4)

axh=fig.add_subplot(gs[1,0])
axh.vlines(inset,0,1,color=NAVY,lw=1.0)
axh.set_xlim(0,N); axh.set_ylim(0,1); axh.set_yticks([]); axh.set_xticks([])
for sp in axh.spines.values(): sp.set_visible(False)
axh.text(-0.012,0.5,"gene set\nmembers",transform=axh.transAxes,ha="right",va="center",fontsize=8.6,color=NAVY)

axr=fig.add_subplot(gs[2,0])
axr.fill_between(rank,0,w,where=w>0,color=RED,alpha=0.5,lw=0)
axr.fill_between(rank,0,w,where=w<=0,color=BLUE,alpha=0.5,lw=0)
axr.set_xlim(0,N); axr.set_ylabel("ranking\nstatistic",fontsize=9)
axr.set_xlabel("all genes, ranked from most up-regulated to most down-regulated")
axr.axhline(0,color="#BBB",lw=1)

# --- ORA ---
ax=fig.add_subplot(gs[:,1])
ax.axis("off")
ax.set_title("ORA — throws most of the data away",loc="left",fontsize=12)
ax.text(0.0,0.90,"① Apply a threshold  (padj < 0.05, |log₂FC| > 1)",fontsize=10.2,color="#333")
ax.text(0.03,0.815,"12,000 genes  →  64 genes",fontsize=11,color=RED,fontweight="bold")
ax.text(0.0,0.73,"② Ask: is my pathway over-represented\n      among those 64, vs. a background set?",fontsize=10.2,
        color="#333",linespacing=1.5,va="top")
ax.text(0.03,0.565,"Fisher's exact test / hypergeometric",fontsize=9.6,color="#555",style="italic")

ax.text(0.0,0.47,"⚠️  Three ways this goes wrong",fontsize=10.6,fontweight="bold",color="#8A3A33")
items=[("Wrong background","using all ~20,000 genes when only\n8,000 were expressed inflates everything"),
       ("Threshold sensitivity","padj < 0.05 vs < 0.10 can change the\ntop pathway entirely"),
       ("Coordinated small changes","50 genes each moving 1.3× — biologically\nreal, invisible to ORA, found by GSEA")]
y=0.37
for t,d in items:
    ax.text(0.03,y,"•  "+t,fontsize=9.8,color="#333",fontweight="bold")
    ax.text(0.065,y-0.030,d,fontsize=9.2,color="#666",linespacing=1.4,va="top")
    y-=0.17
ax.text(0.0,-0.115,"→  Use GSEA when you have a full ranked list.\n     Use ORA when all you have is a gene list from a paper.",
        fontsize=9.8,color="#2E5E4E",linespacing=1.6,va="top")
ax.set_xlim(0,1); ax.set_ylim(-0.16,1)
fig.savefig(OUT+"ora_vs_gsea.png"); print("ok h6")
