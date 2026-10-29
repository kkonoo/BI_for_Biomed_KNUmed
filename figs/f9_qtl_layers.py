from gallery import *
fig,ax=plt.subplots(figsize=(10.5,5.8))
ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis("off")

layers=[("DNA variant","the GWAS hit","",None),
        ("Chromatin accessibility","caQTL","W3",0.10),
        ("DNA methylation","meQTL","W3",0.16),
        ("RNA — expression, splicing","eQTL · sQTL","W4 · W5",0.24),
        ("Protein","pQTL","W6",0.34),
        ("Metabolite","mQTL","",0.46),
        ("Disease","","",1.0)]
x0,x1,h,gap=0.22,0.80,0.095,0.033
for j,(name,qtl,week,a) in enumerate(layers):
    y=0.05+j*(h+gap)
    fc=FAINT if a is None else (CRIMSON if a==1.0 else tint(a))
    ax.add_patch(mp.Rectangle((x0,y),x1-x0,h,fc=fc,ec="none"))
    hot=a==1.0
    ax.text(x0+0.02,y+h/2,name,va="center",fontsize=11.5,fontweight="bold" if hot else "semibold",
            color="white" if hot else TEXT)
    if qtl: ax.text(x1-0.02,y+h/2,qtl,va="center",ha="right",fontsize=11,fontweight="semibold",color=CRIMSON)
    if week: ax.text(x1+0.02,y+h/2,week,va="center",fontsize=10,color=TEXT2)
ytop=0.05+6*(h+gap)
ax.annotate("",xy=(0.13,ytop+h*0.5),xytext=(0.13,0.05+h*0.5),
            arrowprops=dict(arrowstyle="-|>",color="#9AA1AB",lw=1.3))
ax.text(0.11,(ytop+0.05+h)/2,"GWAS\nvariant → disease,\nthe path between\nunseen",ha="right",va="center",fontsize=9.5,color=TEXT2,linespacing=1.4)
ax.text(x1+0.02,ytop+h+0.035,"this course",fontsize=9,color=TEXT2)
header(ax,"Multi-omics QTLs","A variant reaches disease through layers — each layer has its own QTL")
fig.savefig(OUT+"qtl_layers.png"); print("ok f9 qtl")
