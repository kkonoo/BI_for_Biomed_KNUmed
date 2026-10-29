from gallery import *

fig,(ax,bx)=plt.subplots(2,1,figsize=(11,7.0),gridspec_kw={"height_ratios":[1.1,1],"hspace":0.32})

# (a) schematic: a repeat shorter vs. longer than one read
U,F=0.5,2.4                  # repeat-unit width, flank width
def allele(y,n,label):
    ax.add_patch(mp.Rectangle((0,y),F,0.42,fc=FAINT,ec="white",lw=1.2))
    for k in range(n):
        ax.add_patch(mp.Rectangle((F+k*U,y),U,0.42,fc=tint(0.3),ec="white",lw=1.2))
    ax.text(F+n*U/2,y+0.52,f"CAG CAG CAG …  ({n} units)",ha="center",va="bottom",fontsize=9.5,color=CRIMSON)
    ax.add_patch(mp.Rectangle((F+n*U,y),F,0.42,fc=FAINT,ec="white",lw=1.2))
    ax.text(F/2,y+0.21,"unique flank",ha="center",va="center",fontsize=9,color=TEXT2)
    ax.text(F+n*U+F/2,y+0.21,"unique flank",ha="center",va="center",fontsize=9,color=TEXT2)
    ax.text(-0.25,y+0.21,label,ha="right",va="center",fontsize=11,fontweight="semibold",color=TEXT)

RL=6*U+2*1.1                 # one read: spans the 6-unit repeat with 1.1 of flank on each side
def read(x,y,ok,note,mark):
    c=TEXT2 if ok else CRIMSON
    ax.add_patch(mp.FancyArrow(x,y,RL,0,width=0.11,head_width=0.24,head_length=0.2,
                 length_includes_head=True,fc=c,ec="none"))
    ax.text(x+RL+0.35,y+0.11,mark,ha="left",va="center",fontsize=10,color=c,fontweight="bold")
    ax.text(x+RL+0.35,y-0.15,note,ha="left",va="center",fontsize=9.5,color=TEXT2)

allele(2.25,6,"normal allele")
read(F-1.1,1.78,True,"spans the repeat, anchored in both flanks","150 bp read: units can be counted")
allele(0.62,22,"expanded allele")
read(F+3*U,0.15,False,"lies entirely inside the repeat","150 bp read: where is it? how many units?")
ax.set_xlim(-2.4,2*F+22*U+0.3); ax.set_ylim(-0.12,3.0); ax.axis("off")
header(ax,"Tandem repeats","A repeat longer than a read can be neither placed nor counted",x=-0.13)

# (b) real repeat lengths vs. one short read: one hue, darker = closer to disease
rows=[("HTT (CAG)n","Huntington disease",
       [(0,78,tint(0.18),"≤ 26",False),(81,117,tint(0.4),"",False),(120,330,CRIMSON,"≥ 40",True)]),
      ("FMR1 (CGG)n","fragile X",
       [(0,132,tint(0.18),"≤ 44",False),(135,162,tint(0.4),"",False),
        (165,600,tint(0.65),"55–200",False),(603,760,CRIMSON,"> 200",True)])]
for i,(gene,dis,segs) in enumerate(rows):
    y=1-i
    for a,b,c,t,openend in segs:
        bx.add_patch(mp.Rectangle((a,y-0.2),b-a,0.4,fc=c,ec="white",lw=1))
        if openend:
            bx.add_patch(mp.Polygon([(b,y-0.27),(b+22,y),(b,y+0.27)],closed=True,fc=c,ec="none"))
        if t: bx.text((a+b)/2,y,t,ha="center",va="center",fontsize=9.5,fontweight="semibold",
                      color=TEXT if c==tint(0.18) else "white")
    bx.text(-14,y+0.1,gene,ha="right",va="center",fontsize=11,fontweight="semibold",color=TEXT)
    bx.text(-14,y-0.15,dis,ha="right",va="center",fontsize=9.5,color=TEXT2)
bx.text(99,0.66,"27–39",ha="center",va="center",fontsize=8.5,color=TEXT2)
bx.text(148,-0.34,"45–54",ha="center",va="center",fontsize=8.5,color=TEXT2,bbox=dict(fc="white",ec="none",pad=1),zorder=3)
bx.legend(handles=[mp.Patch(fc=c,label=t) for c,t in [(tint(0.18),"normal"),(tint(0.4),"intermediate"),
                                                     (tint(0.65),"premutation (FMR1)"),(CRIMSON,"disease")]],
          title="repeat count",title_fontsize=9,loc="upper right",ncol=4,frameon=False,fontsize=9,
          bbox_to_anchor=(1.0,1.1),handlelength=1.2,columnspacing=1.0,labelcolor=TEXT2)
bx.axvline(150,color=TEXT,ls=(0,(4,3)),lw=1.1)
bx.text(156,1.55,"one short read = 150 bp",ha="left",va="center",fontsize=10,color=TEXT,fontweight="semibold")
bx.set_xlim(0,800); bx.set_ylim(-0.55,1.75)
bx.set_yticks([]); bx.spines["left"].set_visible(False); bx.spines["bottom"].set_color(LINE)
bx.tick_params(colors=TEXT2)
bx.set_xlabel("length of the repeat tract (bp)  =  number of repeats × 3",color=TEXT2)
bx.set_title("Disease alleles are about as long as one read — or far longer",loc="left",
             fontsize=12.5,fontweight="semibold",color=TEXT,x=-0.13,ha="left")
fig.savefig(OUT+"tandem_repeat.png"); print("ok f1 tr")
