from style import *
from scipy.stats import norm
rng=np.random.default_rng(9)
fig,axes=plt.subplots(1,2,figsize=(11.8,4.9),gridspec_kw={"wspace":0.28,"width_ratios":[1.15,1]})

d=0.62
x=np.linspace(-4.2,4.6,600)
ax=axes[0]
ax.fill_between(x,norm.pdf(x,0,1),color=GREY,alpha=0.55,lw=0,label="Controls")
ax.plot(x,norm.pdf(x,0,1),color="#8A8A8A",lw=1.4)
ax.fill_between(x,norm.pdf(x,d,1),color=RED,alpha=0.42,lw=0,label="Cases")
ax.plot(x,norm.pdf(x,d,1),color=RED,lw=1.4)
auc=norm.cdf(d/np.sqrt(2))
ax.set_xlabel("Polygenic risk score  (standardised)"); ax.set_ylabel("density")
ax.set_yticks([])
ax.set_title("A good PRS still overlaps almost completely",loc="left",fontsize=12)
ax.legend(frameon=False,fontsize=10,loc="upper left")
ax.text(0.985,0.95,f"AUC = {auc:.2f}\nmean difference = {d:.2f} SD",transform=ax.transAxes,
        ha="right",va="top",fontsize=10.5,color="#333",linespacing=1.6,
        bbox=dict(boxstyle="round,pad=0.45",fc="#FAFAFA",ec="#DDDDDD"))
ax.annotate("most cases and most controls\nhave indistinguishable scores",
            xy=(0.55,0.235),xytext=(2.75,0.135),fontsize=9.3,color="#444",ha="center",
            arrowprops=dict(arrowstyle="-|>",color="#999",lw=1.2,connectionstyle="arc3,rad=-0.3"))

ax=axes[1]
dec=np.arange(1,11)
z=norm.ppf((dec-0.5)/10)
orr=np.exp(0.62*z); orr=orr/orr[4]*1.0
bars=ax.bar(dec,orr,color=[BLUE]*9+[RED],width=0.72,edgecolor="white")
bars[0].set_color(TEAL)
ax.axhline(1,color="#999",lw=1.1,ls="--")
ax.set_xticks(dec); ax.set_xlabel("PRS decile  (1 = lowest risk)")
ax.set_ylabel("Odds ratio  vs. median decile")
ax.set_title("The same score, shown as deciles",loc="left",fontsize=12)
ax.text(10,orr[-1]+0.06,f"{orr[-1]:.1f}×",ha="center",fontsize=10,fontweight="bold",color=RED)
ax.text(1,orr[0]+0.06,f"{orr[0]:.1f}×",ha="center",fontsize=10,fontweight="bold",color=TEAL)
ax.text(0.5,0.97,"⚠️  Same data, far more impressive.\nExtreme-decile contrasts are the most\ncommon way PRS performance is oversold.",
        transform=ax.transAxes,ha="center",va="top",fontsize=9.3,color="#8A3A33",linespacing=1.5,
        bbox=dict(boxstyle="round,pad=0.45",fc="#FDF3F2",ec="#E8C4C0"))
ax.annotate("",xy=(1,orr[0]+0.14),xytext=(10,orr[-1]+0.30),
            arrowprops=dict(arrowstyle="<|-|>",color="#777",lw=1.2,connectionstyle="arc3,rad=0.18"))
ax.text(5.5,orr[-1]*1.19,f"top vs bottom decile:  {orr[-1]/orr[0]:.1f}×",
        ha="center",fontsize=10,color="#444",style="italic")
ax.set_ylim(0,orr[-1]*1.62)
fig.savefig(OUT+"prs_distribution.png"); print("ok f8")
