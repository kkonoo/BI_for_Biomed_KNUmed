from gallery import *
rng=np.random.default_rng(12)
pos=np.linspace(0,1,90)
fig,axes=plt.subplots(2,2,figsize=(11,5.2),sharex=True,gridspec_kw={"hspace":0.12,"wspace":0.12})

def peak(ax,at,h,hot):
    y=h*np.exp(-((pos-at)/0.08)**2)+rng.uniform(0,0.12,90)
    ax.scatter(pos,y,s=12,color="#C9CED6",lw=0)
    k=np.argmin(np.abs(pos-at))
    if hot: ax.scatter([pos[k]],[y.max()+0.07],s=48,color=CRIMSON,zorder=3)
    else:   ax.scatter([pos[k]],[y.max()+0.07],s=48,facecolor="white",edgecolor=TEXT,lw=1.2,zorder=3)
    ax.set_ylim(0,1.35); ax.set_yticks([])
    for s in ("left","bottom"): ax.spines[s].set_color(LINE)

for col,(gw,eq,hot,title,note) in enumerate([
        (0.45,0.45,True,"Colocalized","same top variant → the GWAS signal\nlikely acts through gene G"),
        (0.38,0.72,False,"Not colocalized","different top variants → two signals that\nhappen to sit in the same locus")]):
    top,bot=axes[0,col],axes[1,col]
    peak(top,gw,1.0,hot); peak(bot,eq,0.85,hot)
    if hot:
        for a in (top,bot): a.axvspan(gw-0.03,gw+0.03,color=FAINT,zorder=0)
    else:
        for a,x in ((top,gw),(bot,eq)): a.axvline(x,color=TEXT2,lw=0.8,ls=(0,(2,2)))
    top.text(0.02,0.93,"GWAS (disease)",transform=top.transAxes,fontsize=9.5,fontweight="semibold",color=TEXT,va="top")
    bot.text(0.02,0.93,"eQTL (gene G)",transform=bot.transAxes,fontsize=9.5,fontweight="semibold",color=TEXT,va="top")
    top.set_title(title,loc="left",fontsize=12,fontweight="bold",color=CRIMSON if hot else TEXT,pad=8)
    bot.set_xticks([]); bot.set_xlabel("position in the locus",color=TEXT2)
    bot.text(0.0,-0.2,note,transform=bot.transAxes,fontsize=9.6,color=TEXT2,va="top",linespacing=1.35)
header(axes[0,0],"Colocalization","Do the disease signal and the expression signal peak at the same variant?",dy=30,x1=2.12)
fig.savefig(OUT+"coloc.png"); print("ok coloc")
