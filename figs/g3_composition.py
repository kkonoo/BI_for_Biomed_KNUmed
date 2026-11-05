from gallery import *
rng=np.random.default_rng(13)
fig,axes=plt.subplots(1,3,figsize=(12.6,4.6),gridspec_kw={"wspace":0.36,"width_ratios":[1,1,1.05]})

ct=["T cell","B cell","Monocyte","NK","Neutrophil"]
# crimson = the cell type that drives the effect; the rest stay greys / light tint
cols=["#5F6673","#9AA1AB","#D5D9DF",tint(0.3),CRIMSON]
ctrl=np.array([0.22,0.09,0.10,0.06,0.53])
case=np.array([0.14,0.06,0.13,0.04,0.63])

# ① composition
ax=axes[0]
bl=np.zeros(2)
for i,(c,col) in enumerate(zip(ct,cols)):
    v=np.array([ctrl[i],case[i]])
    ax.bar([0,1],v,bottom=bl,color=col,width=0.55,edgecolor="white",lw=1.2)
    for k in (0,1):
        if v[k]>=0.08:
            ax.text(k,bl[k]+v[k]/2,c,ha="center",va="center",fontsize=8.8,
                    color="white" if col in (CRIMSON,"#5F6673") else TEXT)
    bl+=v
ax.text(1.33,1-0.63/2,"63%",ha="left",va="center",fontsize=9.5,color=CRIMSON,fontweight="semibold")
ax.text(-0.33,1-0.53/2,"53%",ha="right",va="center",fontsize=9.5,color=CRIMSON,fontweight="semibold")
ax.set_xticks([0,1]); ax.set_xticklabels(["Controls","Cases"])
ax.set_ylabel("cell-type proportion",color=TEXT2); ax.set_ylim(0,1); ax.set_xlim(-0.6,1.6)
ax.text(0.5,-0.2,"colours as in ②",transform=ax.transAxes,ha="center",fontsize=8.5,color=TEXT2)
ax.set_title("① Cases have more neutrophils",loc="left",fontsize=12,fontweight="semibold",color=TEXT,pad=8)

# ② methylation of one CpG by cell type
bx=axes[1]
x=np.arange(5)
meth=np.array([0.22,0.31,0.78,0.28,0.83])
bx.bar(x,meth,color=cols,width=0.62,edgecolor="white",lw=1.2)
bx.set_xticks(x); bx.set_xticklabels(ct,rotation=30,ha="right",fontsize=9.5)
bx.set_ylabel("β value at one CpG",color=TEXT2); bx.set_ylim(0,1)
bx.set_title("② This CpG differs by cell type",loc="left",fontsize=12,fontweight="semibold",color=TEXT,pad=8)

# ③ the bulk tissue: weighted average of ① × ②
cx=axes[2]
sd=0.035
c1=np.dot(ctrl,meth)+rng.normal(0,sd,140)
c2=np.dot(case,meth)+rng.normal(0,sd,140)
parts=cx.violinplot([c1,c2],positions=[0,1],widths=0.62,showmeans=True)
for b,col in zip(parts["bodies"],["#C9CED6",tint(0.45)]):
    b.set_facecolor(col); b.set_alpha(1); b.set_edgecolor("none")
for k in ["cmins","cmaxes","cbars","cmeans"]: parts[k].set_color(TEXT2); parts[k].set_linewidth(1.1)
cx.set_xticks([0,1]); cx.set_xticklabels(["Controls","Cases"])
cx.set_ylabel("β value (whole blood)",color=TEXT2); cx.set_ylim(0.48,0.84)
cx.text(0.5,0.97,"p = 3 × 10⁻⁹",transform=cx.transAxes,ha="center",va="top",fontsize=10.5,
        color=CRIMSON,fontweight="bold")
cx.set_title("③ Apparent differential methylation",loc="left",fontsize=12,fontweight="semibold",color=CRIMSON,pad=8)

for a in axes:
    for s in ("left","bottom"): a.spines[s].set_color(LINE)
    a.tick_params(colors=TEXT2)
header(axes[0],"Cell-type composition · the main confounder",
       "No cell changed its methylation — only the mixture changed",dy=24,x1=3.77)
fig.savefig(OUT+"cell_composition_confounder.png"); print("ok g3")
