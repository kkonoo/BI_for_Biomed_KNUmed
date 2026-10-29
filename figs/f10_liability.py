from gallery import *
from scipy.stats import norm
# K = 5% → threshold T; first-degree relatives of a case: mean shifted by ½·h²·i (h² = 0.6), variance kept at 1 (schematic)
K,h2=0.05,0.6
T=norm.isf(K); i=norm.pdf(T)/K; mu=0.5*h2*i; Kr=norm.sf(T-mu)
x=np.linspace(-3.6,4.6,500)
fig,axes=plt.subplots(2,1,figsize=(10.5,5.6),sharex=True,gridspec_kw={"hspace":0.25})
for ax,m,risk,lab,c in [(axes[0],0,K,"general population",tint(0.55)),(axes[1],mu,Kr,"first-degree relatives of a patient",CRIMSON)]:
    y=norm.pdf(x,m)
    ax.plot(x,y,color=TEXT2,lw=1.3)
    ax.fill_between(x,0,y,where=x<T,color=FAINT)
    ax.fill_between(x,0,y,where=x>=T,color=c)
    ax.text(-3.5,0.33,lab,fontsize=11,fontweight="semibold",color=TEXT,va="center")
    ax.text(T+0.55,0.2,f"affected  {risk*100:.0f}%",fontsize=11,color=c,fontweight="bold",va="center")
    if m:
        ax.axvline(m,color=TEXT2,lw=0.9,ls=(0,(2,2)))
        ax.text(m-0.1,0.5,"shifted right: shares part of the patient's genes",fontsize=9,color=TEXT2,ha="right",va="top")
    ax.set_ylim(0,0.52); ax.set_yticks([])
    for s in ("left",): ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(LINE); ax.tick_params(colors=TEXT2)
for ax in axes:
    ax.axvline(T,color=TEXT,lw=1.3)
axes[0].text(T+0.08,0.47,f"threshold T = {T:.2f} SD  (prevalence K = 5%)",fontsize=9.5,color=TEXT,ha="left",va="center")
axes[1].set_xlabel("liability = genetic + environmental risk, unseen (SD units)",color=TEXT2)
header(axes[0],"Liability-threshold model",
       "A yes/no disease as a continuous liability that crosses a threshold")
fig.savefig(OUT+"liability_threshold.png"); print("ok f10 liab",round(T,2),round(mu,2),round(Kr,3))
