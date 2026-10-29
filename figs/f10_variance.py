from gallery import *
fig,ax=plt.subplots(figsize=(11,4.4))
ax.set_xlim(-0.02,1.02); ax.set_ylim(0,1); ax.axis("off")

# schematic shares of the phenotypic variance V_P
parts=[("V_A",0.50,CRIMSON),("V_D",0.08,tint(0.55)),("V_I",0.05,tint(0.35)),("V_C",0.12,"#C9CED6"),("V_E",0.25,"#E3E6EB")]
x=0; y,h=0.38,0.2
for sym,w,c in parts:
    ax.add_patch(mp.Rectangle((x,y),w,h,fc=c,ec="white",lw=2))
    dark=c in (CRIMSON,tint(0.55))
    ax.text(x+w/2,y+h/2+0.03,f"${sym[0]}_{sym[2]}$",ha="center",va="center",fontsize=14,
            color="white" if dark else TEXT,fontweight="bold")
    x+=w

def bracket(x0,x1,yb,label,c,above=True):
    d=0.04 if above else -0.04
    ax.plot([x0,x0,x1,x1],[yb,yb+d,yb+d,yb],color=c,lw=1.4)
    ax.text(x0,yb+d+(0.03 if above else -0.03),label,ha="left",va="bottom" if above else "top",
            fontsize=11,color=c,fontweight="semibold")
bracket(0,0.50,y+h+0.03,"h² (narrow sense) = $V_A$ / $V_P$   ← what SNPs, GWAS and PRS capture",CRIMSON)
bracket(0,0.63,y+h+0.17,"H² (broad sense) = all genetic $V_G$ / $V_P$",TEXT)
bracket(0,0.63,y-0.03,"genetic  $V_G$",TEXT2,above=False)
bracket(0.63,1.0,y-0.03,"environmental",TEXT2,above=False)
ax.text(0,0.06,"A additive  ·  D dominance  ·  I interaction between genes  ·  C environment shared in a family  ·  E environment unique to each person",
        ha="left",va="center",fontsize=9.5,color=TEXT2)
ax.text(1.0,0.06+0.07,"whole bar = $V_P$, how much the trait varies in this population",ha="right",va="center",fontsize=9.5,color=TEXT2)
header(ax,"Heritability","Split the differences between people into a genetic and an environmental share")
fig.savefig(OUT+"h2_variance.png"); print("ok f10 var")
