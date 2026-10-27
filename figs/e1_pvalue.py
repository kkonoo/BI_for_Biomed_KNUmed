from style import *
from scipy.stats import norm
Z=2.2
P=2*norm.sf(Z)

fig,ax=plt.subplots(figsize=(8.6,3.7))
x=np.linspace(-4.6,4.6,900); y=norm.pdf(x)
ax.plot(x,y,color=NAVY,lw=2.2,zorder=3)
for side in (x>=Z, x<=-Z):
    ax.fill_between(x[side],y[side],color=RED,alpha=0.55,lw=0,zorder=2)
ax.axvline(Z,ymax=0.62,color=RED,lw=2,zorder=4)
ax.axvline(-Z,ymax=0.62,color=RED,lw=1.2,ls="--",alpha=0.6,zorder=4)

ax.text(0,0.43,"distribution of the test statistic\nif H₀ is true",ha="center",fontsize=10.5,
        color=NAVY,fontweight="bold",linespacing=1.4)
ax.annotate(f"observed statistic\nz = {Z}",xy=(Z,0.27),xytext=(3.3,0.33),fontsize=10,color="#8A3A33",
            ha="center",linespacing=1.4,arrowprops=dict(arrowstyle="-|>",color=RED,lw=1.2))
ax.annotate(f"p-value = red area\n= {P:.3f} (both tails)",xy=(2.6,0.008),xytext=(3.55,0.13),fontsize=10,
            color="#8A3A33",ha="center",linespacing=1.4,fontweight="bold",
            arrowprops=dict(arrowstyle="-|>",color=RED,lw=1.2))
ax.annotate("",xy=(-2.6,0.008),xytext=(-3.4,0.11),
            arrowprops=dict(arrowstyle="-|>",color=RED,lw=1.0,alpha=0.6))
ax.text(-3.4,0.12,"just as extreme\nin the other direction",fontsize=9,color="#8A3A33",alpha=0.85,
        ha="center",va="bottom",linespacing=1.35)

ax.set_xlim(-4.6,4.6); ax.set_ylim(0,0.5)
ax.set_xticks([-3,-2,-1,0,1,2,3])
ax.set_xlabel("test statistic (z)")
ax.set_yticks([]); ax.spines["left"].set_visible(False)
fig.savefig(OUT+"pvalue_normal.png"); print("ok e1")
