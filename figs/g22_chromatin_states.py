from gallery import *
from matplotlib.colors import LinearSegmentedColormap

rng=np.random.default_rng(22)
fig,(ax,bx)=plt.subplots(1,2,figsize=(12.4,7.0),gridspec_kw={"width_ratios":[1.0,1.12],"wspace":0.06})

MARKS=["H3K4me3","H3K4me1","H3K36me3","H3K27me3","H3K9me3"]
# Roadmap core 15-state model: (short name, description, RGB, approximate emission per mark)
S=[("TssA","Active TSS",(255,0,0),[.95,.45,.05,.03,.02]),
   ("TssAFlnk","Flanking active TSS",(255,69,0),[.75,.90,.08,.04,.02]),
   ("TxFlnk","Transcription at 5′ / 3′",(50,205,50),[.45,.85,.75,.03,.02]),
   ("Tx","Strong transcription",(0,128,0),[.03,.08,.90,.02,.03]),
   ("TxWk","Weak transcription",(0,100,0),[.02,.10,.35,.02,.03]),
   ("EnhG","Genic enhancer",(194,225,5),[.05,.85,.80,.03,.03]),
   ("Enh","Enhancer",(255,255,0),[.04,.80,.05,.03,.02]),
   ("ZNF/Rpts","ZNF genes & repeats",(102,205,170),[.02,.05,.70,.03,.80]),
   ("Het","Heterochromatin",(138,145,208),[.02,.03,.05,.05,.55]),
   ("TssBiv","Bivalent / poised TSS",(205,92,92),[.85,.50,.04,.80,.05]),
   ("BivFlnk","Flanking bivalent region",(233,150,122),[.40,.85,.05,.75,.05]),
   ("EnhBiv","Bivalent enhancer",(189,183,107),[.05,.65,.04,.70,.04]),
   ("ReprPC","Repressed Polycomb",(128,128,128),[.02,.10,.02,.65,.04]),
   ("ReprPCWk","Weak repressed Polycomb",(192,192,192),[.02,.05,.03,.22,.04]),
   ("Quies","Quiescent / low",(255,255,255),[.01,.02,.02,.02,.02])]
C={n:"#%02X%02X%02X"%rgb for n,_,rgb,_ in S}
cmap=LinearSegmentedColormap.from_list("cr",["#FFFFFF",tint(0.35),CRIMSON])

# left: state list + emission heatmap ------------------------------------------
X0,CW,RH=5.5,0.62,0.5            # heatmap start, cell width, row height
for i,(n,d,rgb,em) in enumerate(S):
    y=-i*RH
    ax.add_patch(mp.Rectangle((0,y-0.16),0.32,0.32,fc=C[n],ec=LINE if n=="Quies" else "none",lw=0.8))
    ax.text(0.45,y,f"{i+1}",ha="left",va="center",fontsize=9.5,color=TEXT2)
    ax.text(0.8,y,n,ha="left",va="center",fontsize=10,fontweight="semibold",color=TEXT)
    ax.text(2.35,y,d,ha="left",va="center",fontsize=9.5,color=TEXT2)
    for j,v in enumerate(em):
        ax.add_patch(mp.Rectangle((X0+j*CW,y-RH/2),CW,RH,fc=cmap(v),ec="white",lw=1.2))
for j,m in enumerate(MARKS):
    ax.text(X0+(j+0.5)*CW,RH*0.62,m,ha="left",va="bottom",rotation=60,fontsize=9.5,
            color=TEXT,fontweight="semibold",rotation_mode="anchor")
for k,(a,b,lab) in enumerate([(0,2,"active TSS"),(2,7,"transcription · enhancers"),
                             (7,9,"heterochromatin"),(9,12,"bivalent"),(12,15,"repressed · quiet")]):
    y0,y1=-a*RH+0.2,-(b-1)*RH-0.2
    if k: ax.plot([0,X0+5*CW],[y0+0.05,y0+0.05],color="#E3E6EB",lw=0.8)
# colour bar
for k,v in enumerate(np.linspace(0,1,40)):
    ax.add_patch(mp.Rectangle((X0+k*(5*CW/40),-15*RH+0.05),5*CW/40,0.17,fc=cmap(v),ec="none"))
ax.text(X0-0.08,-15*RH+0.135,"0",ha="right",va="center",fontsize=8.5,color=TEXT2)
ax.text(X0+5*CW+0.08,-15*RH+0.135,"1",ha="left",va="center",fontsize=8.5,color=TEXT2)
ax.text(X0+2.5*CW,-15*RH-0.2,"how often the mark is present",ha="center",va="top",fontsize=8.5,color=TEXT2)
ax.set_xlim(-0.05,X0+5*CW+0.25); ax.set_ylim(-15*RH-0.6,2.1); ax.axis("off")
header(ax,"Downstream analysis · chromatin states (ChromHMM)",
       "Five noisy tracks become one annotation",dy=26,x1=2.13)
ax.set_title("(a) 15 states, each a combination of marks",loc="left",
             fontsize=12.5,fontweight="semibold",color=TEXT,pad=8)

# right: one locus, five tracks → segmentation ---------------------------------------
x=np.linspace(0,100,1000)
def bump(c,w,h=1.0): return h*np.exp(-0.5*((x-c)/w)**2)
def dom(a,b,h=1.0,e=1.5): return h/(1+np.exp(-(x-a)/e))/(1+np.exp((x-b)/e))
sig={"H3K4me3": bump(22,0.9)+0.06*bump(11,1.2),
     "H3K4me1": bump(11,2.0,0.85)+bump(18.5,1.2,0.7)+bump(25.5,1.3,0.6)+0.35*dom(26,31),
     "H3K36me3": dom(27,55,1.0,1.2)*np.clip((x-20)/32,0.2,1)+0.6*dom(82,86.5,1,0.6),
     "H3K27me3": dom(59,77,0.85,1.2),
     "H3K9me3": dom(82,98,0.8,0.8)}
cols={"H3K4me3":CRIMSON,"H3K4me1":tint(0.6),"H3K36me3":tint(0.4),"H3K27me3":"#7A7F88","H3K9me3":"#4A4E55"}
TH=0.85
for k,m in enumerate(MARKS):
    y0=-k*1.05
    v=np.clip(sig[m]+rng.normal(0,0.035,x.size)+0.04,0,None)*TH
    bx.fill_between(x,y0,y0+v,color=cols[m],lw=0,step="mid")
    bx.plot([0,100],[y0,y0],color=LINE,lw=0.6)
    bx.text(-2,y0+0.3,m,ha="right",va="center",fontsize=10,fontweight="semibold",color=TEXT)
# gene annotation
ya=-5*1.05+0.35
bx.plot([0,100],[ya,ya],color=LINE,lw=1)
bx.add_patch(mp.Rectangle((9,ya-0.13),4,0.26,fc="white",ec=CRIMSON,lw=1))
bx.text(11,ya-0.25,"enhancer",ha="center",va="top",fontsize=9,color=TEXT2)
for a,b in [(22,23.5),(30,31.5),(38,39.5),(46,47.5),(52,55)]:
    bx.add_patch(mp.Rectangle((a,ya-0.16),b-a,0.32,fc=TEXT,ec="none"))
bx.plot([22,55],[ya,ya],color=TEXT,lw=1.6)
bx.annotate("",(25.5,ya+0.5),(22,ya+0.5),arrowprops=dict(arrowstyle="-|>",color=TEXT,lw=1.1))
bx.plot([22,22],[ya+0.16,ya+0.5],color=TEXT,lw=1.1)
bx.text(38.5,ya-0.25,"active gene",ha="center",va="top",fontsize=9,color=TEXT2)
for a,b in [(61,62.5),(67,68.5),(73,75)]:
    bx.add_patch(mp.Rectangle((a,ya-0.16),b-a,0.32,fc="#9A9EA6",ec="none"))
bx.plot([61,75],[ya,ya],color="#9A9EA6",lw=1.6)
bx.text(68,ya-0.25,"silent gene",ha="center",va="top",fontsize=9,color=TEXT2)
bx.add_patch(mp.Rectangle((82,ya-0.13),16,0.26,fc=FAINT,ec=LINE,lw=0.8,hatch="////"))
bx.text(90,ya-0.25,"repeats",ha="center",va="top",fontsize=9,color=TEXT2)
# segmentation bar
SEG=[(0,8,"Quies"),(8,14.5,"Enh"),(14.5,20.5,"TssAFlnk"),(20.5,23.5,"TssA"),(23.5,27,"TssAFlnk"),
     (27,31,"TxFlnk"),(31,51,"Tx"),(51,56,"TxWk"),(56,58.5,"Quies"),(58.5,60,"ReprPCWk"),
     (60,76,"ReprPC"),(76,78.5,"ReprPCWk"),(78.5,82,"Quies"),(82,86.5,"ZNF/Rpts"),(86.5,98,"Het"),(98,100,"Quies")]
ys=ya-1.75
bx.annotate("",(50,ys+0.62),(50,ya-0.62),arrowprops=dict(arrowstyle="-|>",color=TEXT2,lw=1.2))
bx.text(48.5,(ys+0.62+ya-0.62)/2,"ChromHMM",ha="right",va="center",fontsize=9.5,color=TEXT2,style="italic")
for a,b,n in SEG:
    bx.add_patch(mp.Rectangle((a,ys),b-a,0.5,fc=C[n],ec=LINE if n=="Quies" else "white",lw=0.6))
bx.text(-2,ys+0.25,"state",ha="right",va="center",fontsize=10,fontweight="semibold",color=TEXT)
for c,n in [(11.2,"Enh"),(22,"TssA"),(41,"Tx"),(68,"ReprPC"),(92.2,"Het")]:
    bx.text(c,ys-0.12,n,ha="center",va="top",fontsize=9.5,fontweight="semibold",color=TEXT)
bx.text(0,ys-0.75,"schematic, after Roadmap Epigenomics 2015 (core 15-state model, 5 marks)",
        ha="left",va="top",fontsize=8.5,color=TEXT2)
bx.set_xlim(-17,101); bx.set_ylim(ys-1.05,1.25); bx.axis("off")
bx.set_title("(b) One locus, segmented into states",loc="left",fontsize=12.5,
             fontweight="semibold",color=TEXT,pad=8)
fig.savefig(OUT+"chromatin_states.png"); print("ok g22 states")
