from style import *
rng=np.random.default_rng(5)
fig,axes=plt.subplots(1,2,figsize=(11.8,4.8),gridspec_kw={"wspace":0.30,"width_ratios":[1,1.05]})
n=2500
m=rng.normal(0,1.6,n)
p=0.42*m+rng.normal(0,1.45,n)
ax=axes[0]
ax.scatter(m,p,s=6,color=TEAL,alpha=0.32,lw=0)
xx=np.linspace(-5,5,50); ax.plot(xx,0.42*xx,color=NAVY,lw=2.2)
r=np.corrcoef(m,p)[0,1]
ax.set_xlabel("mRNA level  (log scale, centred)"); ax.set_ylabel("protein level  (log scale, centred)")
ax.set_title("Across genes, mRNA explains ~20% of protein",loc="left",fontsize=12)
ax.text(0.965,0.06,f"r = {r:.2f}      r² = {r**2:.2f}",transform=ax.transAxes,va="bottom",ha="right",
        fontsize=12,fontweight="bold",color=NAVY,
        bbox=dict(boxstyle="round,pad=0.4",fc="white",ec="#DDD"))
for x0,y0,t,c in [(-3.4,3.0,"high protein,\nlow mRNA",RED),(3.2,-3.2,"high mRNA,\nlow protein",BLUE)]:
    ax.text(x0,y0,t,fontsize=9.3,color=c,ha="center",linespacing=1.4)
ax.set_xlim(-5,5); ax.set_ylim(-5,5)

ax=axes[1]; ax.axis("off")
ax.set_title("Where the other 80% goes",loc="left",fontsize=12)
rows=[("Translation efficiency","how often a given mRNA is read by ribosomes\n— varies >100-fold between transcripts",TEAL),
      ("Protein half-life","minutes to days. A stable protein accumulates\nfrom a modest transcript",PURPLE),
      ("Post-translational modification","phosphorylation, ubiquitination, cleavage —\nchanges activity with no change in abundance","#C9A227"),
      ("Localisation","the same amount of protein in the nucleus\nversus the cytoplasm is different biology",RED),
      ("Complex stoichiometry","unassembled subunits are degraded, so\nabundance is set by the partner",BLUE)]
y=0.92
for t,d,c in rows:
    ax.plot([0.0,0.022],[y+0.012,y+0.012],color=c,lw=4,solid_capstyle="round")
    ax.text(0.045,y,t,fontsize=10.4,fontweight="bold",color=c,va="top")
    ax.text(0.045,y-0.055,d,fontsize=9.2,color="#555",va="top",linespacing=1.45)
    y-=0.195
ax.set_xlim(0,1); ax.set_ylim(-0.05,1)
fig.text(0.5,-0.03,"None of this is visible to RNA-seq. Whenever a conclusion is about protein, mRNA is a hypothesis — not evidence.",
         ha="center",fontsize=9.6,color="#666")
fig.savefig(OUT+"mrna_vs_protein.png"); print("ok m1")
