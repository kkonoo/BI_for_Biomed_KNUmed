from style import *
rng=np.random.default_rng(37)
fig,axes=plt.subplots(1,2,figsize=(12.4,4.9),gridspec_kw={"wspace":0.30,"width_ratios":[1.25,1]})
methods=["ALDEx2","ANCOM-BC","MaAsLin2","LinDA","DESeq2","edgeR","LEfSe","Wilcoxon\n(CLR)"]
ntax=45
truth=np.zeros(ntax,bool); truth[:8]=True
calls=[]
sens=[0.55,0.70,0.78,0.75,0.88,0.90,0.85,0.62]
fpr =[0.02,0.05,0.10,0.07,0.28,0.34,0.40,0.04]
for s,f in zip(sens,fpr):
    c=np.zeros(ntax,bool)
    c[truth]=rng.random(truth.sum())<s
    c[~truth]=rng.random((~truth).sum())<f
    calls.append(c)
M=np.array(calls)
ax=axes[0]
ax.imshow(M,cmap=plt.matplotlib.colors.ListedColormap(["#F2F2F2",BLUE]),aspect="auto",
          interpolation="nearest")
for j in range(ntax):
    if truth[j]: ax.axvspan(j-0.5,j+0.5,color="#2E8B57",alpha=0.13,lw=0)
ax.set_yticks(range(len(methods))); ax.set_yticklabels(methods,fontsize=9.3)
ax.set_xticks([]); ax.set_xlabel("45 taxa  (green band = the 8 that truly differ)")
ax.set_title("Eight methods, one dataset",loc="left",fontsize=12)
ax.text(1.0,1.02,"blue = called significant",transform=ax.transAxes,fontsize=9,color=BLUE,ha="right")

ax=axes[1]
ncall=M.sum(1); tp=(M&truth).sum(1); fp=(M&~truth).sum(1)
Y=np.arange(len(methods))[::-1]
ax.barh(Y,tp,color="#2E8B57",height=0.62,label="true positives")
ax.barh(Y,fp,left=tp,color=RED,height=0.62,label="false positives",alpha=0.8)
ax.set_yticks(Y); ax.set_yticklabels(methods,fontsize=9.3)
ax.set_xlabel("taxa called significant")
ax.set_title("⚠️ …and they disagree, mostly on false positives",loc="left",fontsize=12,color=RED)
ax.legend(frameon=False,fontsize=9.2,loc="lower right")
agree=(M.sum(0)==len(methods)).sum()
ax.text(0.985,0.97,f"called by all 8 methods:  {agree} taxa\ncalled by exactly one:  {(M.sum(0)==1).sum()} taxa",
        transform=ax.transAxes,ha="right",va="top",fontsize=9.4,color="#333",linespacing=1.5,
        bbox=dict(boxstyle="round,pad=0.42",fc="#FAFAFA",ec="#DDD"))
fig.text(0.5,-0.03,"Report results from more than one method, and treat taxa called by only one as unconfirmed. Simulated here, but this pattern is well documented across dozens of real datasets.",
         ha="center",fontsize=9.4,color="#666")
fig.savefig(OUT+"da_disagreement.png"); print("ok n6")
