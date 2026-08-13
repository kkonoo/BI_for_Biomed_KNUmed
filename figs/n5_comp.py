from style import *
fig,axes=plt.subplots(1,3,figsize=(13.4,4.5),gridspec_kw={"wspace":0.34})
taxa=["Taxon A","Taxon B","Taxon C","Taxon D"]
cols=[BLUE,TEAL,"#C9A227",RED]
abs_c=np.array([100,100,100,20],float)
abs_t=np.array([100,100,100,300],float)

ax=axes[0]
X=np.arange(4); w=0.36
ax.bar(X-w/2,abs_c,w,color=GREY,label="Controls",edgecolor="white")
ax.bar(X+w/2,abs_t,w,color=RED,alpha=0.78,label="Cases",edgecolor="white")
ax.set_xticks(X); ax.set_xticklabels(taxa,fontsize=9)
ax.set_ylabel("absolute abundance (cells / g)")
ax.set_title("① The truth",loc="left",fontsize=11.5)
ax.legend(frameon=False,fontsize=9.2,loc="upper left")
ax.text(0.5,0.80,"only Taxon D changed",transform=ax.transAxes,ha="center",fontsize=10,
        fontweight="bold",color="#2E5E4E")

ax=axes[1]
rc=abs_c/abs_c.sum()*100; rt=abs_t/abs_t.sum()*100
ax.bar(X-w/2,rc,w,color=GREY,edgecolor="white")
ax.bar(X+w/2,rt,w,color=RED,alpha=0.78,edgecolor="white")
ax.set_xticks(X); ax.set_xticklabels(taxa,fontsize=9)
ax.set_ylabel("relative abundance (%)")
ax.set_title("② ⚠️ What sequencing measures",loc="left",fontsize=11.5,color=RED)
for i in range(4):
    fc=rt[i]/rc[i]
    ax.text(i,max(rc[i],rt[i])+1.6,f"{fc:.2f}×",ha="center",fontsize=9.2,fontweight="bold",
            color=RED if abs(np.log2(fc))>0.2 else "#888")
ax.set_ylim(0,max(rc.max(),rt.max())*1.30)
ax.text(0.5,0.03,"A, B and C all appear to 'decrease'\nbecause D took up more of the total.",
        transform=ax.transAxes,ha="center",fontsize=9.3,color="#8A3A33",linespacing=1.5,
        bbox=dict(boxstyle="round,pad=0.4",fc="#FDF3F2",ec="#E8C4C0"))

ax=axes[2]
def clr(x):
    lx=np.log(x); return lx-lx.mean()
cc,ct=clr(rc),clr(rt)
ax.bar(X-w/2,cc,w,color=GREY,edgecolor="white")
ax.bar(X+w/2,ct,w,color=TEAL,alpha=0.85,edgecolor="white")
ax.axhline(0,color="#999",lw=1.1)
ax.set_xticks(X); ax.set_xticklabels(taxa,fontsize=9)
ax.set_ylabel("CLR-transformed abundance")
ax.set_title("③ After CLR transformation",loc="left",fontsize=11.5,color=TEAL)
ax.text(0.5,0.99,"Ratios between taxa within a sample are\npreserved. CLR-based tests (ALDEx2,\nANCOM-BC) work on these ratios rather\nthan on raw proportions.",
        transform=ax.transAxes,ha="center",va="top",fontsize=9.2,color="#2E5E4E",linespacing=1.5,
        bbox=dict(boxstyle="round,pad=0.4",fc="#F1F7F4",ec="#BFD8CC"))
ax.set_ylim(-1.65,1.35)
fig.suptitle("⚠️  Compositional data: sequencing measures proportions, and proportions must sum to 1",
             x=0.075,ha="left",fontsize=13,fontweight="bold",y=1.04)
fig.savefig(OUT+"compositional.png"); print("ok n5")
