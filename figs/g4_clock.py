from style import *
rng=np.random.default_rng(17)
fig,axes=plt.subplots(1,2,figsize=(11.5,4.8),gridspec_kw={"wspace":0.28,"width_ratios":[1.1,1]})
n=260
age=rng.uniform(20,85,n)
pred=age+rng.normal(0,4.0,n)
ax=axes[0]
ax.scatter(age,pred,s=17,color=TEAL,alpha=0.55,lw=0)
lim=[15,92]; ax.plot(lim,lim,color="#999",ls="--",lw=1.3)
ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel("Chronological age (years)"); ax.set_ylabel("Predicted (epigenetic) age")
ax.set_title("An epigenetic clock",loc="left",fontsize=12)
ax.text(0.04,0.95,"r = 0.96\nmedian error ≈ 3–4 years",transform=ax.transAxes,va="top",
        fontsize=10,color="#333",linespacing=1.5)
i=np.argmax(pred-age); j=np.argmin(pred-age)
for k,c,lab in [(i,RED,"age acceleration\n(predicted > actual)"),(j,BLUE,"age deceleration")]:
    ax.plot([age[k],age[k]],[age[k],pred[k]],color=c,lw=2.0,zorder=3)
    ax.scatter([age[k]],[pred[k]],s=60,color=c,zorder=4,ec="white",lw=1.2)
    ax.annotate(lab,xy=(age[k],pred[k]),xytext=(0,18 if k==i else -34),
                textcoords="offset points",ha="center",fontsize=9,color=c,linespacing=1.3)
ax.text(0.985,0.05,"the residual — not the prediction —\nis what people study",transform=ax.transAxes,
        ha="right",fontsize=9.3,color="#555",style="italic",linespacing=1.4)

ax=axes[1]
gens=["1st gen\n(Horvath, Hannum)","2nd gen\n(PhenoAge, GrimAge)","3rd gen\n(DunedinPACE)"]
trained=["chronological\nage","age + clinical\nbiomarkers","rate of change\nover time"]
pred_of=[0.30,0.68,0.72]
bars=ax.barh(range(3),pred_of,color=[GREY,OLIVE,TEAL],height=0.52,edgecolor="white")
ax.set_yticks(range(3)); ax.set_yticklabels(gens,fontsize=9.5)
ax.invert_yaxis(); ax.set_xlim(0,1.0); ax.set_xticks([])
ax.set_xlabel("relative association with mortality / morbidity")
ax.set_title("Not all clocks measure the same thing",loc="left",fontsize=12)
for k,t in enumerate(trained):
    ax.text(pred_of[k]+0.03,k,"trained on: "+t.replace("\n"," "),va="center",fontsize=8.8,color="#555")
ax.text(0.5,-0.30,"⚠️  A clock trained to predict chronological age is, by construction,\npenalised for detecting biological ageing.",
        transform=ax.transAxes,ha="center",va="top",fontsize=9.4,color="#8A3A33",linespacing=1.5,
        bbox=dict(boxstyle="round,pad=0.45",fc="#FDF3F2",ec="#E8C4C0"))
fig.savefig(OUT+"epigenetic_clock.png"); print("ok g4")
