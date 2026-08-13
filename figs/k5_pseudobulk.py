from style import *
from scipy.stats import norm
rng=np.random.default_rng(66)
fig,axes=plt.subplots(1,3,figsize=(13.4,4.5),gridspec_kw={"wspace":0.32,"width_ratios":[1.1,1,1]})

ND=6
donor_eff=rng.normal(0,0.55,ND)
grp=np.array([0,0,0,1,1,1])
true=0.10
ax=axes[0]
allv=[];cols=[]
for d in range(ND):
    n=900
    v=rng.normal(3.0+donor_eff[d]+true*grp[d],0.75,n)
    allv.append(v); cols.append(GREY if grp[d]==0 else RED)
    ax.scatter(np.full(n,d)+rng.normal(0,0.10,n),v,s=3,color=cols[-1],alpha=0.28,lw=0)
    ax.plot([d-0.30,d+0.30],[v.mean()]*2,color="#222",lw=2.4,zorder=5)
ax.set_xticks(range(ND)); ax.set_xticklabels([f"C{i+1}" for i in range(3)]+[f"T{i+1}" for i in range(3)])
ax.set_ylabel("expression of GENE-Y per cell")
ax.set_title("5,400 cells from 6 donors",loc="left",fontsize=12)
ax.text(0.02,0.97,"black bars = donor means\n(they differ a lot)",transform=ax.transAxes,va="top",
        fontsize=9.3,color="#333",linespacing=1.4)

ax=axes[1]
cells=np.concatenate(allv); g=np.repeat(grp,900)
from scipy.stats import ttest_ind
t1=ttest_ind(cells[g==1],cells[g==0])
pb=np.array([v.mean() for v in allv])
t2=ttest_ind(pb[grp==1],pb[grp==0])
labs=["per-cell t-test\n(n = 5,400)","pseudobulk\n(n = 6 donors)"]
ps=[t1.pvalue,t2.pvalue]
b=ax.bar([0,1],[-np.log10(p) for p in ps],color=[RED,TEAL],width=0.5,edgecolor="white")
ax.axhline(-np.log10(0.05),color="#888",ls="--",lw=1.2)
ax.text(-0.42,-np.log10(0.05)+3.0,"p = 0.05",ha="left",fontsize=9,color="#777")
ax.set_xticks([0,1]); ax.set_xticklabels(labs,fontsize=9.5)
ax.set_ylabel(r"$-\log_{10}(p)$")
ax.set_title("Same data, two tests",loc="left",fontsize=12)
for i,p in enumerate(ps):
    ax.text(i,-np.log10(p)+4.5,f"p = {p:.0e}" if p<0.01 else f"p = {p:.2f}",
            ha="center",fontsize=10.5,fontweight="bold",color=[RED,TEAL][i])
ax.set_ylim(0,max(-np.log10(ps[0])*1.28,4))

ax=axes[2]
ax.axis("off")
ax.set_title("Why the left-hand test is wrong",loc="left",fontsize=12)
ax.text(0,0.86,"Cells from one donor are NOT\nindependent observations.",fontsize=10.6,color="#333",
        linespacing=1.5,va="top")
ax.text(0,0.66,"A per-cell test treats 900 cells from\none person as 900 independent\nsamples. It is measuring how many\ncells you loaded, not how many\npeople you studied.",
        fontsize=10,color="#444",linespacing=1.55,va="top")
ax.text(0,0.29,"→  The unit of replication is the\n     DONOR, not the cell.",
        fontsize=10.8,color="#2E5E4E",fontweight="bold",linespacing=1.55,va="top")
ax.text(0,0.10,"Sum counts within each donor × cell type,\nthen use DESeq2 / edgeR exactly as in W4.",
        fontsize=9.6,color="#555",style="italic",linespacing=1.5,va="top")
ax.set_xlim(0,1); ax.set_ylim(0,1)
fig.suptitle("⚠️  Pseudobulk: the single most important statistical point in single-cell analysis",
             x=0.065,ha="left",fontsize=13,fontweight="bold",y=1.03)
fig.savefig(OUT+"pseudobulk.png"); print("ok k5")
