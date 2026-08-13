from style import *
rng=np.random.default_rng(41)
fig,axes=plt.subplots(1,2,figsize=(11.8,4.7),gridspec_kw={"wspace":0.28})
P=4000
inten=rng.normal(20,2.6,P)
pmiss=1/(1+np.exp((inten-18.2)*1.5))
miss=pmiss*100+rng.normal(0,4,P); miss=np.clip(miss,0,100)
ax=axes[0]
ax.scatter(inten,miss,s=6,color=PURPLE,alpha=0.30,lw=0)
o=np.argsort(inten)
ax.plot(inten[o],pmiss[o]*100,color=NAVY,lw=2.4)
ax.set_xlabel("mean log₂ intensity (when detected)"); ax.set_ylabel("% of runs where the protein is missing")
ax.set_title("Missingness is not random",loc="left",fontsize=12)
ax.axvspan(12,17.5,color=RED,alpha=0.10,lw=0)
ax.text(15.0,62,"MNAR\nbelow the limit\nof detection",ha="center",fontsize=9.4,color="#8A3A33",linespacing=1.45)
ax.axvspan(21.5,29,color=TEAL,alpha=0.10,lw=0)
ax.text(24.2,40,"MAR / MCAR\nstochastic DDA\nselection",ha="center",fontsize=9.4,color="#2E5E4E",linespacing=1.45)
ax.set_xlim(12,29); ax.set_ylim(-4,104)

ax=axes[1]
n=1200
lfc=rng.normal(0,0.30,n); se=rng.uniform(0.16,0.35,n)
z=lfc/se; p=2*np.exp(-0.5*z**2)*0.5
nm=260
lfc_i=np.concatenate([lfc,rng.normal(2.7,0.45,nm)])
p_i=np.concatenate([p,10**(-rng.uniform(2.2,8.5,nm))])
ax.scatter(lfc,-np.log10(p),s=7,color="#C8C8C8",alpha=0.6,lw=0,label="genuinely measured")
ax.scatter(lfc_i[n:],-np.log10(p_i[n:]),s=10,color=RED,alpha=0.7,lw=0,
           label="created by imputing zeros")
ax.axhline(-np.log10(0.05),color="#888",ls="--",lw=1.1)
ax.set_xlabel("log₂ fold change"); ax.set_ylabel("−log₁₀(p)")
ax.set_title("⚠️ What naive imputation produces",loc="left",fontsize=12,color=RED)
ax.legend(frameon=False,fontsize=9.2,loc="upper left")
ax.text(0.985,0.06,"Replacing missing values with a small constant\nturns 'not detected in group A' into a large,\nhighly significant fold change — for every\nprotein near the detection limit.",
        transform=ax.transAxes,ha="right",fontsize=9.1,color="#333",linespacing=1.5,
        bbox=dict(boxstyle="round,pad=0.42",fc="#FDF3F2",ec="#E8C4C0"))
fig.savefig(OUT+"missing_values.png"); print("ok m6")
