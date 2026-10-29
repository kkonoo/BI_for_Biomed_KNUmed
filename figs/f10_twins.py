from gallery import *
rng=np.random.default_rng(4)
fig=plt.figure(figsize=(12.5,4.8))
gs=fig.add_gridspec(1,3,width_ratios=[1,1.05,1.1],wspace=0.35,top=0.8)
a1,a2,a3=[fig.add_subplot(gs[0,k]) for k in range(3)]

def sub(ax,num,title,q):
    x0=ax.get_position().x0; top=a1.get_position().y1
    fig.text(x0,top+0.085,f"{num}  {title}",fontsize=12,fontweight="bold",color=TEXT,va="bottom")
    fig.text(x0,top+0.03,q,fontsize=9.3,color=TEXT2,va="bottom")

# 1. genes shared with a relative
rel=[("MZ twin",1.0),("parent, child,\nsibling, DZ twin",0.5),("grandparent · uncle/aunt",0.25),("first cousin",0.125)]
for j,(name,f) in enumerate(rel):
    y=len(rel)-1-j
    a1.barh(y,1,color=FAINT,height=0.55); a1.barh(y,f,color=CRIMSON if j==0 else tint(0.35+0.5*f),height=0.55)
    a1.text(-0.04,y,name,ha="right",va="center",fontsize=9.6,color=TEXT,linespacing=1.15)
    a1.text(f+0.03 if f<0.8 else f-0.03,y,f"{f*100:g}%",ha="left" if f<0.8 else "right",va="center",fontsize=10,
            fontweight="semibold",color=TEXT if f<0.8 else "white")
a1.set_xlim(0,1.05); a1.axis("off")
sub(a1,"1","Genes shared","the closer the relative, the more alleles in common")

# 2. twin pairs: MZ resemble each other more than DZ
for r,c,lab,z in [(0.45,"#9AA1AB","DZ pairs  r = 0.45",1),(0.80,CRIMSON,"MZ pairs  r = 0.80",2)]:
    xy=rng.multivariate_normal([0,0],[[1,r],[r,1]],140)
    a2.scatter(xy[:,0],xy[:,1],s=12,color=c,alpha=0.75,lw=0,zorder=z,label=lab)
a2.legend(frameon=False,fontsize=9.5,loc="upper left",handletextpad=0.2,markerscale=1.6,labelcolor=TEXT)
a2.set_xlabel("twin 1 (trait, SD)",color=TEXT2); a2.set_ylabel("twin 2",color=TEXT2)
a2.set_xlim(-3.2,3.2); a2.set_ylim(-3.2,3.2); a2.set_xticks([]); a2.set_yticks([])
for s in ("left","bottom"): a2.spines[s].set_color(LINE)
sub(a2,"2","Twin correlation","same home; only the share of genes differs")

# 3. the arithmetic (Falconer)
a3.set_xlim(0,1); a3.set_ylim(0,1); a3.axis("off")
rows=[("MZ:  r = A + C","0.80",[(0.70,CRIMSON),(0.10,"#C9CED6")]),
      ("DZ:  r = ½A + C","0.45",[(0.35,tint(0.5)),(0.10,"#C9CED6")])]
for j,(eq,val,segs) in enumerate(rows):
    y=0.78-j*0.24; x=0.0
    a3.text(0,y+0.075,eq,fontsize=11,color=TEXT,fontweight="semibold",va="bottom")
    for w,c in segs:
        a3.add_patch(mp.Rectangle((x,y),w*0.8,0.065,fc=c,ec="white",lw=1.5)); x+=w*0.8
    a3.text(x+0.02,y+0.032,val,va="center",fontsize=10.5,color=TEXT)
a3.text(0,0.24,"difference = ½A = 0.35",fontsize=10.5,color=TEXT2,va="center")
a3.text(0,0.1,"h² = 2 × (0.80 − 0.45) = 0.70",fontsize=13,color=CRIMSON,fontweight="bold",va="center")
for k,(c,t) in enumerate([(CRIMSON,"A  genetic"),("#C9CED6","C  shared environment")]):
    a3.add_patch(mp.Rectangle((k*0.3,0.36),0.03,0.03,fc=c,ec="none"))
    a3.text(k*0.3+0.045,0.375,t,fontsize=9,color=TEXT2,va="center")
sub(a3,"3","Falconer's formula","twice the extra MZ resemblance = heritability")

fig.text(0.125,1.08,"Estimating heritability from relatives",fontsize=10.5,color=CRIMSON)
fig.text(0.125,1.025,"MZ twins resemble each other more than DZ twins — that excess is genetic",fontsize=16,fontweight="bold",color=TEXT)
fig.savefig(OUT+"twin_design.png"); print("ok f10 twins")
