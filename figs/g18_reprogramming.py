from gallery import *

fig,ax=plt.subplots(figsize=(11.5,6.0))
DARK="#4A4F57"; IMP="#8A8F98"
def curve(pts,n=400,sm=9):
    xs,ys=zip(*pts); x=np.linspace(xs[0],xs[-1],n); y=np.interp(x,xs,ys)
    k=np.ones(sm)/sm; yp=np.pad(y,(sm//2,sm//2),mode="edge"); return x,np.convolve(yp,k,mode="valid")

# wave bands
for a,b,lab,sub in [(0.0,3.7,"Wave 1 — early embryo","genome-wide loss; imprinted regions are protected"),
                    (5.8,8.6,"Wave 2 — germ line (PGCs)","near-complete erasure, imprints too; then reset by sex")]:
    ax.add_patch(mp.Rectangle((a,0),b-a,1.0,fc=FAINT,ec="none",zorder=0))
    ax.text(a+0.12,1.07,lab,ha="left",va="bottom",fontsize=11.5,fontweight="bold",color=CRIMSON)
    ax.text(a+0.12,1.02,sub,ha="left",va="bottom",fontsize=9.5,color=TEXT2)

# wave 1: paternal fast (active), maternal slow (passive), then de novo re-methylation
xp,yp=curve([(0,0.86),(0.35,0.80),(0.9,0.36),(1.6,0.30),(2.6,0.21),(3.2,0.22),(3.9,0.55),(4.6,0.74)],sm=15)
xm,ym=curve([(0,0.42),(0.9,0.40),(1.8,0.30),(2.6,0.21),(3.2,0.22),(3.9,0.55),(4.6,0.74)],sm=15)
ax.plot(xp[xp<=2.7],yp[xp<=2.7],color=CRIMSON,lw=2.4,zorder=3)
ax.plot(xm[xm<=2.7],ym[xm<=2.7],color=tint(0.45),lw=2.4,zorder=3)
xs,ys=curve([(2.6,0.21),(3.2,0.22),(3.9,0.55),(4.6,0.74),(5.9,0.74),(6.7,0.55),(7.4,0.2),(8.0,0.07),(8.5,0.06)],sm=15)
ax.plot(xs,ys,color=DARK,lw=2.4,zorder=3)
# gametogenesis: sex-specific re-methylation
xs1,ys1=curve([(8.4,0.06),(8.8,0.2),(9.4,0.78),(10.0,0.85),(11.0,0.86)],sm=15)
xs2,ys2=curve([(8.4,0.06),(9.6,0.07),(10.2,0.2),(10.7,0.4),(11.0,0.42)],sm=15)
ax.plot(xs1,ys1,color=CRIMSON,lw=2.4,zorder=3); ax.plot(xs2,ys2,color=tint(0.45),lw=2.4,zorder=3)
# imprinted regions: one allele methylated (~50%); protected in wave 1, erased + reset in wave 2
xi,yi=curve([(0,0.5),(6.0,0.5),(6.8,0.36),(7.6,0.08),(8.5,0.05),(9.3,0.12),(10.2,0.5),(11.0,0.5)],sm=15)
ax.plot(xi,yi,color=IMP,lw=1.8,ls=(0,(4,3)),zorder=2)

# annotations
kw=dict(fontsize=9.5,ha="left",va="center")
ax.text(0.95,0.7,"paternal genome:\nfast, active loss (TET3)",color=CRIMSON,**kw)
ax.text(1.25,0.13,"maternal genome: slow, passive loss\nover cell divisions",color=TEXT,**kw)
ax.annotate("de novo re-methylation\n(DNMT3A/B)",(3.85,0.53),(4.15,0.3),fontsize=9.5,color=TEXT,ha="left",va="center",
            arrowprops=dict(arrowstyle="-",color=TEXT2,lw=0.9))
ax.text(1.0,0.535,"imprinted regions stay ~50%",color=IMP,fontsize=9.2,ha="left",va="bottom",fontweight="semibold")
ax.annotate("near-complete\nerasure (~5–10%)",(8.0,0.07),(6.55,0.17),fontsize=9.5,color=TEXT,ha="right",va="center",
            arrowprops=dict(arrowstyle="-",color=TEXT2,lw=0.9))
ax.text(9.55,0.92,"sperm: re-methylated before birth",color=CRIMSON,fontsize=9.5,ha="left",va="center")
ax.text(10.5,0.17,"oocyte: during\noocyte growth",color=TEXT,fontsize=9.5,ha="left",va="center")
ax.text(10.75,0.6,"imprints reset\nby sex",color=IMP,fontsize=9,ha="left",va="center",fontweight="semibold")

# stage axis
stages=[(0,"gametes /\nfertilisation"),(0.9,"zygote"),(2.45,"blastocyst"),(3.6,"implantation"),
        (4.85,"epiblast /\nsomatic"),(6.0,"PGCs\nspecified"),(7.15,"migrating\nPGCs"),(8.2,"gonad"),
        (9.6,"gameto-\ngenesis"),(11.0,"gametes")]
for x,t in stages:
    ax.plot([x,x],[0,-0.02],color=TEXT2,lw=0.9,clip_on=False)
    ax.text(x,-0.04,t,ha="center",va="top",fontsize=9.2,color=TEXT,linespacing=1.15)
ax.annotate("",(11.95,-0.055),(11.55,-0.055),arrowprops=dict(arrowstyle="-|>",color=TEXT2,lw=1))
ax.text(12.0,-0.055,"next\ngeneration",ha="left",va="center",fontsize=9,color=TEXT2)

ax.set_xlim(-0.35,12.8); ax.set_ylim(-0.16,1.15)
ax.set_yticks([0,0.25,0.5,0.75,1]); ax.set_yticklabels(["0%","25%","50%","75%","100%"])
ax.tick_params(colors=TEXT2,labelsize=9); ax.spines["left"].set_bounds(0,1); ax.spines["left"].set_color(LINE)
ax.spines["bottom"].set_visible(False); ax.set_xticks([])
ax.set_ylabel("genome-wide CpG methylation",color=TEXT2,fontsize=10)
# legend
lg=[(CRIMSON,"-","paternal / male germ line"),(tint(0.45),"-","maternal / female germ line"),
    (DARK,"-","embryo and early germ line"),(IMP,(0,(4,3)),"imprinted regions")]
xl=[0.0,2.9,6.0,8.9]
for k,(c,ls,t) in enumerate(lg):
    ax.plot([xl[k],xl[k]+0.4],[-0.235,-0.235],color=c,lw=2.2 if k<3 else 1.8,ls=ls,clip_on=False)
    ax.text(xl[k]+0.5,-0.235,t,ha="left",va="center",fontsize=9,color=TEXT2)
header(ax,"Epigenetic reprogramming in mammals","DNA methylation is wiped twice per generation",dy=8)
ax.text(-0.35,-0.31,"A mark acquired by a parent has to escape both waves → the reason transgenerational inheritance "
        "is contested in mammals.",ha="left",va="top",fontsize=10.5,color=TEXT,fontweight="semibold")
ax.text(-0.35,-0.365,"Schematic based on mouse data; the human timeline is similar but slower.",ha="left",va="top",
        fontsize=9,color=TEXT2)
fig.savefig(OUT+"reprogramming.png"); print("ok g18 reprogramming")
