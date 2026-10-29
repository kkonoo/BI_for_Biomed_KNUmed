from gallery import *
rng=np.random.default_rng(7)
fig,ax=plt.subplots(figsize=(12.5,5.6))
ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis("off")

W=np.array([0.55,0,0.3,0,-0.25,0.1])          # learned weights for 6 SNPs near gene G
GREY={0:"#F2F3F5",1:"#C9CED6",2:"#8A919C"}
cw,ch,nr=0.026,0.058,7                          # cell width/height, people shown

def grid(x0,y0,G):
    for i in range(G.shape[0]):
        for j in range(G.shape[1]):
            ax.add_patch(mp.Rectangle((x0+j*cw,y0-i*ch),cw*0.92,ch*0.88,fc=GREY[G[i,j]],ec="none"))
    for j in range(G.shape[1]):
        ax.text(x0+j*cw+cw*0.46,y0+ch*1.15,str(j+1),ha="center",fontsize=8,color=TEXT2)
    ax.text(x0,y0+ch*2.0,"SNP",fontsize=8.5,color=TEXT2)

def column(x0,y0,vals,label,dashed=False):
    lo,hi=vals.min(),vals.max()
    for i,v in enumerate(vals):
        ax.add_patch(mp.Rectangle((x0,y0-i*ch),cw*1.1,ch*0.88,fc=tint(0.15+0.75*(v-lo)/(hi-lo)),
                     ec=TEXT2 if dashed else "none",lw=0.8,ls=(0,(2,1.5)) if dashed else "-"))
    ax.text(x0+cw*0.55,y0+ch*1.15,label,ha="center",fontsize=8.5,color=CRIMSON,linespacing=1.15,va="bottom")

def weights(x0,y0,alpha=1):
    ax.plot([x0,x0+6*cw],[y0,y0],color=LINE,lw=0.8)
    for j,w in enumerate(W):
        if w: ax.add_patch(mp.Rectangle((x0+j*cw+cw*0.15,y0),cw*0.62,w*0.12,fc=TEXT,ec="none",alpha=alpha))

def step(x,num,name,sub,hot=False):
    ax.scatter([x+0.015],[0.9],s=420,color=CRIMSON if hot else tint(0.16),zorder=3)
    ax.text(x+0.015,0.9,num,ha="center",va="center",fontsize=10,fontweight="bold",color="white" if hot else CRIMSON,zorder=4)
    ax.text(x+0.045,0.915,name,fontsize=12,fontweight="bold",color=TEXT,va="center")
    ax.text(x+0.045,0.885,sub,fontsize=9,color=TEXT2,va="top",linespacing=1.3)

# 1. learn: reference panel with genotypes and measured expression
step(0.0,"1","Learn","reference set (e.g. GTEx):\ngenotypes + measured expression")
G1=rng.integers(0,3,(nr,6)); e1=G1@W+rng.normal(0,0.15,nr)
grid(0.02,0.63,G1); column(0.2,0.63,e1,"measured\nexpression of G")
weights(0.02,0.16)
ax.text(0.02,0.095,"learned weight per SNP",fontsize=9,color=TEXT2)
ax.text(0.02,0.05,"expression ≈ Σ weight × genotype",fontsize=9.5,color=TEXT,fontweight="semibold")

# 2. predict: GWAS cohort with genotypes only
step(0.35,"2","Predict","GWAS cohort:\ngenotypes only, no expression")
G2=rng.integers(0,3,(nr,6)); e2=G2@W
grid(0.37,0.63,G2); column(0.55,0.63,e2,"predicted\nexpression of G",dashed=True)
weights(0.37,0.16,alpha=0.45)
ax.text(0.37,0.095,"same weights",fontsize=9,color=TEXT2)
ax.text(0.37,0.05,"→ every person gets a predicted value",fontsize=9.5,color=TEXT,fontweight="semibold")

# 3. test: predicted expression vs disease, gene by gene
step(0.69,"3","Test","predicted expression\nvs. disease",hot=True)
ia=fig.add_axes([0.75,0.3,0.17,0.42])
ctrl=rng.normal(0,1,120); case=rng.normal(0.55,1,120)
ia.scatter(rng.uniform(-0.17,0.17,120),ctrl,s=9,color="#C9CED6",lw=0)
ia.scatter(1+rng.uniform(-0.17,0.17,120),case,s=9,color=tint(0.6),lw=0)
for x,v,c in ((0,ctrl,TEXT2),(1,case,CRIMSON)): ia.plot([x-0.25,x+0.25],[v.mean()]*2,color=c,lw=2.4)
ia.set_xticks([0,1],["controls","cases"]); ia.set_yticks([]); ia.set_xlim(-0.55,1.55)
ia.set_ylabel("predicted expression of G",color=TEXT2,fontsize=9)
for s in ("left","bottom"): ia.spines[s].set_color(LINE)
ia.tick_params(colors=TEXT2,labelsize=9.5)
ax.text(0.71,0.095,"repeat for every gene (~20,000 tests,",fontsize=9,color=TEXT2)
ax.text(0.71,0.05,"not millions of variants)",fontsize=9,color=TEXT2)
header(ax,"TWAS — transcriptome-wide association study","Predict each person's gene expression from genotype, then test it against the disease")
fig.savefig(OUT+"twas.png"); print("ok twas")
