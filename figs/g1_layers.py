from gallery import *
W,H=11.0,5.0; Y0,Y1=0.0,1.0
ASP=((Y1-Y0)/H)/(1.0/W)          # y-units per x-unit for a visually round circle
fig,ax=plt.subplots(figsize=(W,H))
def circ(x,y,r,**kw): ax.add_patch(mp.Ellipse((x,y),2*r,2*r*ASP,**kw))
def square(x,y,r,**kw): ax.add_patch(mp.Rectangle((x-r,y-r*ASP),2*r,2*r*ASP,**kw))

REP="#5F6673"; K4=CRIMSON; K27=tint(0.45)       # repressive mark · H3K4me3 · H3K27ac
NUC_FC="#E3E6EB"; NUC_EC="#9AA1AB"; DNA="#9AA1AB"

def mark(x,y,kind,r=0.0105):
    if kind=="rep": square(x,y,r*0.9,fc=REP,ec="white",lw=0.8,zorder=4)
    else:           circ(x,y,r,fc=K4 if kind=="k4" else K27,ec="white",lw=0.8,zorder=4)

def nucleosome(x,y,kinds,r=0.030):
    circ(x,y,r,fc=NUC_FC,ec=NUC_EC,lw=1.1,zorder=3)
    for dx,k in zip((-0.0145,0.0145),kinds):
        y1=y+r*ASP*0.55; y2=y1+0.055
        ax.plot([x,x+dx],[y1,y2],color=NUC_EC,lw=1.0,zorder=2)
        mark(x+dx,y2,k)

y0=0.70
# closed side
ax.text(0.175,0.965,"Closed chromatin",ha="center",fontsize=12.5,fontweight="bold",color=TEXT)
ax.text(0.175,0.915,"gene silenced",ha="center",fontsize=9.5,color=TEXT2)
ax.plot([0.03,0.32],[y0,y0],color=DNA,lw=2.0,zorder=1)
for x in np.linspace(0.055,0.295,5): nucleosome(x,y0,("rep","rep"))
ax.text(0.175,y0-0.105,"nucleosomes packed tightly",ha="center",fontsize=9.5,color=TEXT2)

# open side
ax.text(0.70,0.965,"Open chromatin",ha="center",fontsize=12.5,fontweight="bold",color=CRIMSON)
ax.text(0.70,0.915,"gene active",ha="center",fontsize=9.5,color=CRIMSON)
ax.plot([0.46,0.99],[y0,y0],color=DNA,lw=2.0,zorder=1)
for x in [0.492,0.556,0.836,0.900,0.964]: nucleosome(x,y0,("k4","k27"))
ax.add_patch(mp.FancyBboxPatch((0.625,y0-0.062),0.155,0.124,
             boxstyle="round,pad=0.004,rounding_size=0.015",fc="white",ec=CRIMSON,lw=1.3,ls=(0,(3,2)),zorder=1))
ax.text(0.7025,y0+0.115,"nucleosome-free region",ha="center",fontsize=9.5,color=TEXT)
ax.annotate("",xy=(0.7025,y0-0.078),xytext=(0.7025,y0-0.145),
            arrowprops=dict(arrowstyle="-|>",color=CRIMSON,lw=1.4))
ax.text(0.7025,y0-0.20,"Tn5 inserts here → ATAC-seq peak",ha="center",fontsize=9.5,color=CRIMSON,
        fontweight="semibold")

# DNA methylation on the same sequence
yl=0.36
ax.plot([0.03,0.97],[yl,yl],color=LINE,lw=1.6,zorder=1)
for xs,filled in [(np.linspace(0.055,0.295,8),True),(np.linspace(0.635,0.770,5),False)]:
    for x in xs:
        ax.plot([x,x],[yl,yl+0.048],color=NUC_EC,lw=1.1,zorder=2)
        circ(x,yl+0.055,0.0092,fc=TEXT if filled else "white",ec=TEXT,lw=1.3,zorder=3)
ax.text(0.175,yl-0.055,"DNA methylation high  →  promoter off",ha="center",fontsize=9.5,color=TEXT)
ax.text(0.7025,yl-0.055,"DNA methylation low  →  promoter on",ha="center",fontsize=9.5,color=TEXT)

# legend
ax.plot([0.03,0.97],[0.235,0.235],color="#E3E6EB",lw=0.9)
leg=[("rep","H3K27me3 / H3K9me3 — repressive",0.055,0.165),
     ("k4","H3K4me3 — active promoter",0.055,0.110),
     ("k27","H3K27ac — active promoter or enhancer",0.055,0.055)]
for k,t,x,y in leg:
    mark(x,y,k); ax.text(x+0.025,y,t,va="center",fontsize=9.5,color=TEXT)
circ(0.60,0.165,0.0092,fc=TEXT,ec=TEXT,lw=1.3)
ax.text(0.625,0.165,"methylated cytosine",va="center",fontsize=9.5,color=TEXT)
circ(0.60,0.110,0.0092,fc="white",ec=TEXT,lw=1.3)
ax.text(0.625,0.110,"unmethylated cytosine",va="center",fontsize=9.5,color=TEXT)

ax.set_xlim(0,1); ax.set_ylim(Y0,Y1); ax.axis("off")
header(ax,"The epigenome","Same DNA sequence, different state — everything that differs here is epigenetic",dy=4)
fig.savefig(OUT+"epigenome_layers.png"); print("ok g1")
