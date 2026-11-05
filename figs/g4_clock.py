from gallery import *
rng=np.random.default_rng(17)
fig,axes=plt.subplots(1,2,figsize=(11.5,4.9),gridspec_kw={"wspace":0.42,"width_ratios":[1.0,1.0]})

# (a) predicted vs chronological age; the residual is the quantity of interest
n=260
age=rng.uniform(20,85,n)
pred=age+rng.normal(0,4.0,n)
ax=axes[0]
lim=[15,92]; ax.plot(lim,lim,color=TEXT2,ls=(0,(4,3)),lw=1.1)
ax.scatter(age,pred,s=15,color="#C9CED6",lw=0)
ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel("chronological age (years)",color=TEXT2); ax.set_ylabel("predicted (epigenetic) age",color=TEXT2)
ax.text(0.97,0.04,"r = 0.96 · median error ≈ 3–4 years\nresidual = predicted − actual\n→ the quantity clock studies analyse",
        transform=ax.transAxes,ha="right",va="bottom",fontsize=9.5,color=TEXT2,linespacing=1.45)
i=np.argmax(pred-age); j=np.argmin(pred-age)
for k,c,lab,dx,ha in [(i,CRIMSON,"age acceleration\npredicted > actual",-10,"right"),
                      (j,TEXT,"age deceleration\npredicted < actual",10,"left")]:
    ax.plot([age[k],age[k]],[age[k],pred[k]],color=c,lw=2.0,zorder=3)
    ax.scatter([age[k]],[pred[k]],s=55,color=c,zorder=4,ec="white",lw=1.2)
    ax.annotate(lab,xy=(age[k],pred[k]),xytext=(dx,0),textcoords="offset points",ha=ha,va="center",
                fontsize=9.5,color=c,fontweight="semibold",linespacing=1.3,
                bbox=dict(fc="white",ec="none",pad=1))
ax.set_title("An epigenetic clock",loc="left",fontsize=12.5,fontweight="semibold",color=TEXT,pad=8)

# (b) clock generations: what they were trained on → how well they track health
bx=axes[1]
gens=["1st gen\nHorvath, Hannum","2nd gen\nPhenoAge, GrimAge","3rd gen\nDunedinPACE"]
trained=["chronological age","age + clinical biomarkers,\nmortality","rate of ageing\nover time"]
bx.text(0.33,0.30,"→ penalised, by design,\n   for tracking biological ageing",va="top",fontsize=9,color=CRIMSON,linespacing=1.3)
val=[0.30,0.68,0.72]
cols=["#C9CED6",tint(0.45),CRIMSON]
bx.barh(range(3),val,color=cols,height=0.5,edgecolor="white")
bx.set_yticks(range(3)); bx.set_yticklabels(gens,fontsize=9.5,linespacing=1.3)
bx.invert_yaxis(); bx.set_xlim(0,1.0); bx.set_xticks([])
for k,t in enumerate(trained):
    bx.text(val[k]+0.03,k,"trained on\n"+t,va="center",fontsize=9,color=TEXT2,linespacing=1.3)
bx.set_xlabel("relative association with mortality / morbidity",color=TEXT2)
bx.set_title("Not all clocks measure the same thing",loc="left",fontsize=12.5,fontweight="semibold",color=TEXT,pad=8)

for a in axes:
    for s in ("left","bottom"): a.spines[s].set_color(LINE)
    a.tick_params(colors=TEXT2)
bx.tick_params(axis="y",colors=TEXT)
header(axes[0],"Epigenetic clocks","The residual, not the predicted age, carries the biology",dy=24,x1=2.42)
fig.savefig(OUT+"epigenetic_clock.png"); print("ok g4")
