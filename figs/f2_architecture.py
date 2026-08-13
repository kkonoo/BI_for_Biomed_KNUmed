from style import *
fig,ax=plt.subplots(figsize=(9,6))
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(3e-4,0.62); ax.set_ylim(0.93,90)
ax.set_xlabel("Allele frequency  (log scale)")
ax.set_ylabel("Effect size,  odds ratio  (log scale)")

f=np.logspace(np.log10(3e-4),np.log10(0.62),400)
for n,c,lab,lx in [(2e3,"#C4C4C4","detection limit, n = 2,000",3.0e-3),
                   (5e5,"#909090","detection limit, n = 500,000",3.0e-3)]:
    ax.plot(f,1+14/np.sqrt(n*f),color=c,lw=2.0,zorder=1)
    ax.text(lx,1+14/np.sqrt(n*lx),"  "+lab,fontsize=9,color="#7A7A7A",va="bottom",ha="left")

pts=[(6e-4,25,RED,"Mendelian\ndisease",1.35e-3,52),
     (2.2e-2,2.8,PURPLE,"Low-frequency\nvariants",4.0e-2,7.5),
     (0.22,1.12,BLUE,"Common complex\ndisease",0.10,1.75)]
for x,y,c,t,tx,ty in pts:
    ax.scatter([x],[y],s=340,color=c,alpha=0.25,zorder=3)
    ax.scatter([x],[y],s=70,color=c,zorder=4)
    ax.annotate(t,xy=(x,y),xytext=(tx,ty),fontsize=11,fontweight="bold",color=c,
                ha="center",va="bottom",zorder=5,
                arrowprops=dict(arrowstyle="-",color=c,lw=1.1,alpha=0.55,
                                shrinkA=2,shrinkB=8))

ax.text(3.6e-4,0.945,"rare AND small effect — essentially undetectable at any realistic n",
        fontsize=9.5,color="#AAAAAA",ha="left",va="bottom")

ax.text(0.985,0.965,
        "How each is found\n"
        "•  Mendelian → family linkage, rare-variant sequencing\n"
        "•  Low-frequency → imputation + large biobanks\n"
        "•  Common → GWAS (§3)",
        transform=ax.transAxes,fontsize=9,color="#444",ha="right",va="top",
        multialignment="left", linespacing=1.6,
        bbox=dict(boxstyle="round,pad=0.5",fc="#FAFAFA",ec="#DDDDDD"))

ax.set_xticks([1e-3,1e-2,1e-1,0.5]); ax.set_xticklabels(["0.1%","1%","10%","50%"])
ax.set_yticks([1,2,5,10,50]); ax.set_yticklabels(["1.0","2","5","10","50"])
ax.set_title("Allelic architecture: which study design finds which variants",loc="left")
fig.text(0.5,-0.02,"Grey curves are approximate. Anything plotting below a curve cannot be detected at that sample size.",
         ha="center",fontsize=9,color="#666")
fig.savefig(OUT+"allelic_architecture.png"); print("ok")
