from gallery import *
from matplotlib.patches import Arc
rng=np.random.default_rng(5)
# chr16, GRCh37 (Mb). Gene bodies simplified; causal variant rs1421085 in FTO intron 1
genes=[("RPGRIP1L",53.633,53.738,-1),("FTO",53.738,54.155,1),("IRX3",54.317,54.322,-1),("IRX5",54.965,54.970,1)]
snp=53.801
fig,ax=plt.subplots(figsize=(11.5,5.6))
ax.set_xlim(53.55,55.08); ax.set_ylim(-1.62,2.35); ax.axis("off")

# a small association tower around the variant
p=np.sort(rng.uniform(53.70,53.92,160)); d=np.abs(p-snp)
y=0.55+1.2*np.exp(-(d/0.025)**2)*rng.uniform(0.7,1,160)+rng.uniform(0,0.12,160)
ax.scatter(p,y,s=10,color="#C9CED6",lw=0)
ax.scatter([snp],[y.max()+0.08],s=60,color=CRIMSON,zorder=3)
ax.text(snp+0.014,y.max()+0.08,"rs1421085",ha="left",va="center",fontsize=10.5,color=CRIMSON,fontweight="semibold")
ax.text(snp,y.max()+0.35,"obesity GWAS signal",ha="center",va="center",fontsize=9.5,color=TEXT2)

# gene track
ax.plot([53.55,55.08],[0,0],color=LINE,lw=1.2)
for name,s,e,strand in genes:
    w=max(e-s,0.012)
    ax.add_patch(mp.Rectangle((s,-0.13),w,0.26,fc=FAINT,ec=TEXT2,lw=0.9,zorder=2))
    lab={"RPGRIP1L":(s+w/2,-0.3,"top"),"FTO":(54.0,0.2,"bottom")}.get(name,(s+w/2,0.2,"bottom"))
    ax.text(lab[0],lab[1],name,ha="center",va=lab[2],fontsize=11,fontstyle="italic",fontweight="semibold",color=TEXT)
for x in [53.738,53.811,53.89,53.97,54.03,54.07,54.10,54.13,54.153]:   # FTO exons (schematic)
    ax.add_patch(mp.Rectangle((x,-0.13),0.004,0.26,fc=TEXT2,ec="none",zorder=3))
ax.plot([snp,snp],[0.16,0.5],color=CRIMSON,lw=1,ls=(0,(2,2)))

# long-range enhancer contacts (below the track)
for tgt in (54.3195,54.9675):
    c=(snp+tgt)/2; wd=tgt-snp
    ax.add_patch(Arc((c,-0.15),wd,1.6*wd,theta1=180,theta2=360,color=CRIMSON,lw=1.6))
ax.text(54.06,-0.62,"0.5 Mb",ha="center",va="top",fontsize=9.5,color=CRIMSON)
ax.text(54.384,-1.14,"1.2 Mb",ha="center",va="top",fontsize=9.5,color=CRIMSON)

# the two errors, as numbered notes
def note(x,y,num,head,body):
    ax.scatter([x],[y],s=380,color=tint(0.16),zorder=3)
    ax.text(x,y,num,ha="center",va="center",fontsize=10,fontweight="bold",color=CRIMSON,zorder=4)
    ax.text(x+0.028,y+0.09,head,ha="left",va="center",fontsize=10.5,fontweight="semibold",color=TEXT)
    ax.text(x+0.028,y-0.15,body,ha="left",va="center",fontsize=9.5,color=TEXT2)
note(54.42,1.75,"1","“The hit is in the FTO gene”","it sits in an intron — non-coding, so it acts by regulation")
note(54.42,1.05,"2","“So FTO is the causal gene”","the enhancer switches on IRX3 and IRX5, 0.5–1.2 Mb away")
ax.text(54.31,-1.5,"in fat-cell precursors: risk allele → more IRX3/IRX5 → fat stored instead of burned (Claussnitzer et al. 2015)",
        ha="center",fontsize=9,color=TEXT2)
header(ax,"Example: the FTO obesity locus","The hit lies in an FTO intron — the genes it acts on are IRX3 and IRX5")
fig.savefig(OUT+"fto_locus.png"); print("ok f9 fto")
