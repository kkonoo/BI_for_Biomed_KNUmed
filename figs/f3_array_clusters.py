from gallery import *
rng=np.random.default_rng(3)

def cloud(n,deg,sd=2.6):
    th=np.deg2rad(rng.normal(deg,sd,n)); r=rng.normal(1,0.11,n)
    return r*np.cos(th), r*np.sin(th)

fig,axes=plt.subplots(1,2,figsize=(11,5.0),gridspec_kw={"wspace":0.18})
cols={"AA":"#D5D9DF","AB":"#9AA1AB","BB":"#5F6673"}   # greys: crimson is kept for the rare carrier

# (a) common SNP: three clean clusters (MAF 30%, HWE: 49 / 42 / 9 %)
ax=axes[0]
for g,deg,n in [("AA",5,490),("AB",45,420),("BB",85,90)]:
    x,y=cloud(n,deg); ax.scatter(x,y,s=11,color=cols[g],lw=0,alpha=0.9)
    lx,ly={"AA":(0.66,0.08),"AB":(0.93,0.93),"BB":(0.36,1.1)}[g]
    ax.text(lx,ly,g,ha="right" if g=="AA" else "center",va="center",fontsize=12,fontweight="bold",color=TEXT)
x,y=cloud(10,np.array([24]*6+[66]*4),sd=4)
ax.scatter(x,y,s=16,facecolor="white",edgecolor="#9AA1AB",lw=0.9)
ax.text(0.98,0.98,"common SNP (MAF 30%)\n1,000 people",transform=ax.transAxes,ha="right",va="top",fontsize=10,color=TEXT2)
ax.annotate("between clusters → no call",xy=(x[0],y[0]),xytext=(1.0,0.6),fontsize=9.5,color=TEXT2,
            arrowprops=dict(arrowstyle="-",color="#9AA1AB",lw=0.8))

# (b) very rare variant: one carrier among 1,000, plus two noisy non-carriers
bx=axes[1]
x,y=cloud(997,5); bx.scatter(x,y,s=11,color=cols["AA"],lw=0)
bx.text(0.66,0.08,"AA",ha="right",va="center",fontsize=12,fontweight="bold",color=TEXT)
bx.scatter([0.62],[0.6],s=70,color=CRIMSON,zorder=3)
bx.annotate("the one real carrier (AB)\n— no cluster to compare with",xy=(0.62,0.6),xytext=(0.12,1.2),
            fontsize=9.5,color=CRIMSON,fontweight="semibold",arrowprops=dict(arrowstyle="-",color=CRIMSON,lw=0.8))
nx,ny=[0.93,1.02],[0.36,0.28]
bx.scatter(nx,ny,s=30,facecolor="white",edgecolor=TEXT2,lw=1.1,zorder=3)
bx.annotate("noisy AA samples:\ncalled as carriers?",xy=(nx[0],ny[0]),xytext=(1.0,0.72),fontsize=9.5,color=TEXT2,
            arrowprops=dict(arrowstyle="-",color="#9AA1AB",lw=0.8))
bx.text(0.98,0.98,"very rare variant\n1 carrier in 1,000 people",transform=bx.transAxes,ha="right",va="top",fontsize=10,color=TEXT2)

for a in axes:
    a.set_xlim(-0.08,1.5); a.set_ylim(-0.08,1.5); a.set_aspect("equal")
    a.set_xticks([]); a.set_yticks([])
    for s in ("left","bottom"): a.spines[s].set_color(LINE)
    a.set_xlabel("allele A signal",color=TEXT2,fontsize=10); a.set_ylabel("allele B signal",color=TEXT2,fontsize=10)
header(axes[0],"Array genotype calling","Genotypes are read off clusters — a very rare variant has none")
fig.savefig(OUT+"array_clusters.png"); print("ok f3 clusters")
