from style import *
import matplotlib.patches as mp
W,H=11.0,5.2; Y0,Y1=0.05,1.0
ASP=((Y1-Y0)/H)/(1.0/W)          # y-units per x-unit for a visually round circle
fig,ax=plt.subplots(figsize=(W,H))
def circ(x,y,r,**kw): ax.add_patch(mp.Ellipse((x,y),2*r,2*r*ASP,**kw))

def nucleosome(x,y,r=0.030,marks=None):
    circ(x,y,r,fc="#B9C4CF",ec="#7A8894",lw=1.2,zorder=3)
    for dx,dy,c in (marks or []):
        ax.plot([x,x+dx],[y+r*ASP*0.55,y+r*ASP*0.55+dy],color="#7A8894",lw=1.0,zorder=2)
        circ(x+dx,y+r*ASP*0.55+dy,0.0105,fc=c,ec="white",lw=0.8,zorder=4)

y0=0.70
ax.text(0.175,0.955,"Closed chromatin",ha="center",fontsize=12.5,fontweight="bold",color="#6B7280")
ax.text(0.175,0.905,"gene silenced",ha="center",fontsize=9.5,color="#6B7280",style="italic")
ax.plot([0.03,0.32],[y0,y0],color="#9AA5B1",lw=2.2,zorder=1)
for x in np.linspace(0.055,0.295,5):
    nucleosome(x,y0,marks=[(-0.020,0.055,"#7B4F9E"),(0.020,0.055,"#7B4F9E")])
ax.text(0.175,y0-0.105,"nucleosomes packed tightly",ha="center",fontsize=9,color="#666")

ax.text(0.70,0.955,"Open chromatin",ha="center",fontsize=12.5,fontweight="bold",color=TEAL)
ax.text(0.70,0.905,"gene active",ha="center",fontsize=9.5,color=TEAL,style="italic")
ax.plot([0.46,0.97],[y0,y0],color="#9AA5B1",lw=2.2,zorder=1)
for x in [0.492,0.556,0.845,0.909,0.958]:
    nucleosome(x,y0,marks=[(-0.020,0.055,"#2E8B57"),(0.020,0.055,"#C9A227")])
ax.add_patch(mp.FancyBboxPatch((0.625,y0-0.062),0.155,0.124,
             boxstyle="round,pad=0.004,rounding_size=0.015",fc="#E7F3F3",ec=TEAL,lw=1.6,ls="--",zorder=1))
ax.text(0.7025,y0+0.115,"nucleosome-free region",ha="center",fontsize=9,color=TEAL)
ax.annotate("",xy=(0.7025,y0-0.078),xytext=(0.7025,y0-0.145),
            arrowprops=dict(arrowstyle="-|>",color=TEAL,lw=1.5))
ax.text(0.7025,y0-0.205,"Tn5 inserts here → ATAC-seq peak",ha="center",fontsize=9,color=TEAL)

yl=0.36
ax.plot([0.03,0.97],[yl,yl],color="#C9C9C9",lw=1.6,zorder=1)
for xs,filled in [(np.linspace(0.055,0.295,8),True),(np.linspace(0.635,0.770,5),False)]:
    for x in xs:
        ax.plot([x,x],[yl,yl+0.048],color="#8A8A8A",lw=1.1,zorder=2)
        circ(x,yl+0.055,0.0092,fc="#C1544B" if filled else "white",ec="#C1544B",lw=1.4,zorder=3)
ax.text(0.175,yl-0.055,"DNA methylation high  →  promoter off",ha="center",fontsize=9.3,color="#C1544B")
ax.text(0.7025,yl-0.055,"DNA methylation low  →  promoter on",ha="center",fontsize=9.3,color="#C1544B")

leg=[("H3K27me3 / H3K9me3 — repressive","#7B4F9E",0.055,0.185),
     ("H3K4me3 — active promoter","#2E8B57",0.055,0.130),
     ("H3K27ac — active enhancer","#C9A227",0.055,0.075)]
for t,c,x,y in leg:
    circ(x,y,0.0105,fc=c,ec="white",lw=0.8); ax.text(x+0.025,y,t,va="center",fontsize=9.5,color="#333")
circ(0.60,0.185,0.0092,fc="#C1544B",ec="#C1544B",lw=1.4)
ax.text(0.625,0.185,"methylated cytosine",va="center",fontsize=9.5,color="#333")
circ(0.60,0.130,0.0092,fc="white",ec="#C1544B",lw=1.4)
ax.text(0.625,0.130,"unmethylated cytosine",va="center",fontsize=9.5,color="#333")

ax.set_xlim(0,1); ax.set_ylim(Y0,Y1); ax.axis("off")
ax.set_title("The epigenome: same sequence, different state",loc="left",x=0.0,y=1.01)
fig.text(0.5,0.0,"The DNA sequence is identical on both sides. Everything that differs here is epigenetic.",
         ha="center",fontsize=9.5,color="#666")
fig.savefig(OUT+"epigenome_layers.png"); print("ok")
