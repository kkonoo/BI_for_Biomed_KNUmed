from gallery import *
# height: share of phenotypic variance explained
rows=[("Twin studies","all genetic effects",0.80,None,"#C9CED6"),
      ("All common SNPs together","SNP heritability",0.45,None,tint(0.35)),
      ("Significant SNPs, 2010","180 loci, ~180,000 people",0.10,None,tint(0.55)),
      ("Significant SNPs, 2022","12,111 SNPs, 5.4 million people — European ancestry",0.40,None,CRIMSON),
      ("    same SNPs, other ancestries","",0.10,0.20,tint(0.25))]
fig,ax=plt.subplots(figsize=(11,4.6))
for j,(name,sub,v,v2,c) in enumerate(rows):
    y=len(rows)-1-j
    ax.barh(y,(v+v2)/2 if v2 else v,color=c if not v2 else tint(0.45),height=0.56)
    if v2: ax.errorbar((v+v2)/2,y,xerr=(v2-v)/2,color=TEXT2,lw=1.1,capsize=4)
    ax.text(-0.015,y+(0.12 if sub else 0),name,ha="right",va="center",fontsize=10.5,color=TEXT,fontweight="semibold")
    if sub: ax.text(-0.015,y-0.2,sub,ha="right",va="center",fontsize=8.8,color=TEXT2)
    lab=f"{v*100:.0f}–{v2*100:.0f}%" if v2 else f"~{v*100:.0f}%"
    ax.text((v2 or v)+0.012,y,lab,va="center",fontsize=10.5,color=CRIMSON if c==CRIMSON else TEXT,fontweight="semibold")
ax.annotate("",xy=(0.10,2.38),xytext=(0.45,2.38),arrowprops=dict(arrowstyle="<->",color=TEXT2,lw=1))
ax.text(0.275,2.48,"the \"missing\" part, 2010",ha="center",va="bottom",fontsize=9,color=TEXT2)
ax.set_xlim(0,1); ax.set_ylim(-0.6,len(rows)-0.4)
ax.set_yticks([]); ax.spines["left"].set_visible(False); ax.spines["bottom"].set_color(LINE)
ax.set_xticks([0,0.25,0.5,0.75,1],["0","25%","50%","75%","100%"]); ax.tick_params(colors=TEXT2)
ax.set_xlabel("share of height variance explained",color=TEXT2)
header(ax,"Missing heritability — height","Bigger GWAS closed most of the gap — in European ancestry",x=-0.34)
fig.savefig(OUT+"missing_heritability.png"); print("ok f10 missing")
