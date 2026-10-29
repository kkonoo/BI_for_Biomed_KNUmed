from gallery import *
rng=np.random.default_rng(9)
fig=plt.figure(figsize=(12.5,4.9))
gs=fig.add_gridspec(2,3,height_ratios=[1,1],width_ratios=[1,1.15,1.25],wspace=0.32,hspace=0.12,top=0.8)
a1=fig.add_subplot(gs[:,0]); a2=fig.add_subplot(gs[0,1]); a3=fig.add_subplot(gs[1,1],sharex=a2)
a4=fig.add_subplot(gs[:,2])

def sub(ax,num,title,q):
    """Panel title in figure coords, so all three line up whatever the panel height."""
    x0=ax.get_position().x0; top=a1.get_position().y1
    fig.text(x0,top+0.085,f"{num}  {title}",fontsize=12,fontweight="bold",color=TEXT,va="bottom")
    fig.text(x0,top+0.03,q,fontsize=9.3,color=TEXT2,va="bottom")

# 1. eQTL: expression by genotype at the GWAS variant
for g,(lab,m) in enumerate([("AA",3.2),("AB",5.6),("BB",8.1)]):
    v=rng.normal(m,1.1,60)
    a1.scatter(g+rng.uniform(-0.18,0.18,60),v,s=12,color="#C9CED6",lw=0)
    a1.plot([g-0.25,g+0.25],[m,m],color=CRIMSON,lw=2.4)
a1.plot([0,1,2],[3.2,5.6,8.1],color=CRIMSON,lw=1,ls=(0,(3,2)))
a1.set_xticks([0,1,2],["AA","AB","BB"]); a1.set_yticks([])
a1.set_xlabel("genotype at the GWAS variant",color=TEXT2); a1.set_ylabel("expression of gene G",color=TEXT2)
sub(a1,"1","eQTL","Does the variant change a gene's expression?")

# 2. colocalization: GWAS and eQTL peaks share one variant
pos=np.linspace(0,1,90); peak=0.52
for ax,lab,h in [(a2,"GWAS (disease)",1.0),(a3,"eQTL (gene G)",0.85)]:
    y=h*np.exp(-((pos-peak)/0.09)**2)+rng.uniform(0,0.12,90)
    ax.scatter(pos,y,s=12,color="#C9CED6",lw=0)
    k=np.argmin(np.abs(pos-peak)); ax.scatter([pos[k]],[y.max()+0.07],s=46,color=CRIMSON,zorder=3)
    ax.axvspan(peak-0.03,peak+0.03,color=FAINT,zorder=0)
    ax.text(0.02,0.93,lab,transform=ax.transAxes,fontsize=9.5,color=TEXT,va="top",fontweight="semibold")
    ax.set_yticks([]); ax.set_ylim(0,1.35)
plt.setp(a2.get_xticklabels(),visible=False); a3.set_xticks([])
a3.set_xlabel("position in the locus",color=TEXT2)
sub(a2,"2","Colocalization","Do the two signals share one causal variant?")

# 3. TWAS in three steps
a4.set_xlim(0,1); a4.set_ylim(0,1); a4.axis("off")
steps=[("Learn","in a reference set with genotypes and expression (e.g. GTEx):","which SNPs predict gene G's expression?"),
       ("Predict","in the GWAS cohort (genotypes only):","each person's genetically predicted expression"),
       ("Test","predicted expression ~ disease,","one test per gene (~20,000, not millions)")]
for j,(name,l1,l2) in enumerate(steps):
    y=0.8-j*0.33
    a4.scatter([0.05],[y],s=420,color=CRIMSON if j==2 else tint(0.16),zorder=3)
    a4.text(0.05,y,str(j+1),ha="center",va="center",fontsize=10,fontweight="bold",color="white" if j==2 else CRIMSON,zorder=4)
    a4.text(0.13,y+0.055,name,fontsize=11,fontweight="semibold",color=TEXT,va="center")
    a4.text(0.13,y-0.005,l1,fontsize=9.2,color=TEXT2,va="center")
    a4.text(0.13,y-0.06,l2,fontsize=9.2,color=TEXT2,va="center")
a4.plot([0.05,0.05],[0.14,0.8],color=LINE,lw=1.1,ls=(0,(2,2.5)),zorder=1)
sub(a4,"3","TWAS","Does genetically predicted expression track the disease?")

for ax in (a1,a2,a3):
    for s in ("left","bottom"): ax.spines[s].set_color(LINE)
    ax.tick_params(colors=TEXT2)
fig.text(0.125,1.08,"From variant to gene",fontsize=10.5,color=CRIMSON)
fig.text(0.125,1.025,"Three questions that point a GWAS hit at its gene",fontsize=16,fontweight="bold",color=TEXT)
fig.savefig(OUT+"variant_to_gene.png"); print("ok f9 v2g")
