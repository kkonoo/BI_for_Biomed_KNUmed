from gallery import *

fig,ax=plt.subplots(figsize=(11.5,6.6))
MONO="DejaVu Sans Mono"; GREY="#8A8F98"

def base(x,y,ch,fc,tc="white",w=0.36,h=0.34,ec="none"):
    ax.add_patch(mp.FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle="round,pad=0.0,rounding_size=0.06",
                                   fc=fc,ec=ec,lw=0.9,zorder=4))
    ax.text(x,y,ch,ha="center",va="center",fontsize=11,family=MONO,fontweight="bold",color=tc,zorder=5)

def scenario(y,title,sub,target,ext,col,channel):
    ax.text(0.0,y+1.05,title,ha="left",va="center",fontsize=11.5,fontweight="semibold",color=TEXT)
    ax.text(0.0,y+0.76,sub,ha="left",va="center",fontsize=9,color=TEXT2)
    ax.add_patch(mp.Circle((0.55,y-0.1),0.42,fc=FAINT,ec=LINE,lw=1,zorder=2))
    ax.text(0.55,y-0.1,"bead",ha="center",va="center",fontsize=8.5,color=TEXT2,zorder=3)
    ax.plot([0.97,4.95],[y-0.1,y-0.1],color=TEXT2,lw=3,solid_capstyle="butt",zorder=3)    # probe
    ax.text(2.9,y-0.4,"50-mer probe",ha="center",va="center",fontsize=8.8,color=TEXT2)
    ax.text(4.95,y-0.38,"3′",ha="center",va="center",fontsize=8.5,color=TEXT2)
    ax.plot([1.15,6.4],[y+0.32,y+0.32],color=TEXT,lw=1.6,zorder=3)                        # target DNA
    ax.text(6.5,y+0.32,"bisulfite-converted DNA",ha="left",va="center",fontsize=8.8,color=TEXT2)
    for xx in np.arange(1.25,4.9,0.22): ax.plot([xx,xx],[y-0.06,y+0.28],color=LINE,lw=0.8,zorder=1)
    base(5.25,y+0.32,target,col if target=="C" else "white",tc="white" if target=="C" else TEXT,
         ec="none" if target=="C" else GREY)
    ax.text(5.25,y+0.62,"target",ha="center",va="bottom",fontsize=8.5,color=TEXT2)
    base(5.25,y-0.1,ext,col if ext=="G" else GREY)
    ax.plot([5.25,5.25],[y+0.08,y+0.15],color=LINE,lw=0.8)
    ax.add_patch(mp.RegularPolygon((5.25,y-0.52),5,radius=0.13,fc=col if ext=="G" else GREY,ec="none",zorder=4))
    ax.annotate(f"labelled {ext} added → {channel} channel",(5.42,y-0.52),(5.75,y-0.52),fontsize=9.8,
                fontweight="semibold",color=col if ext=="G" else TEXT,ha="left",va="center",
                arrowprops=dict(arrowstyle="-",color=LINE,lw=0.9))

ax.text(0.0,5.25,"① Infinium II: one probe per CpG, one extra base",ha="left",va="center",fontsize=12,
        fontweight="semibold",color=TEXT)
scenario(3.75,"Methylated CpG","5mC was protected → target is still C","C","G",CRIMSON,"green")
scenario(1.55,"Unmethylated CpG","C became U → T after bisulfite","T","A",GREY,"red")

# ② intensities → β
RX=9.35
ax.plot([RX-0.35,RX-0.35],[0.25,5.45],color="#E3E6EB",lw=0.9)
ax.text(RX,5.25,"② Intensities → β",ha="left",va="center",fontsize=12,fontweight="semibold",color=TEXT)
ax.text(RX,4.85,"each bead carries many probe copies;\neach copy reports one DNA molecule",ha="left",va="center",
        fontsize=8.8,color=TEXT2,linespacing=1.4)
for k,(lab,v,c) in enumerate([("M  (methylated)",0.8,CRIMSON),("U  (unmethylated)",0.2,GREY)]):
    y=4.05-k*0.62
    ax.text(RX,y+0.2,lab,ha="left",va="center",fontsize=9.5,color=TEXT)
    ax.add_patch(mp.Rectangle((RX,y-0.17),2.2*v,0.2,fc=c,ec="none"))
ax.text(RX,2.45,"β = M / (M + U + 100)",ha="left",va="center",fontsize=12.5,fontweight="bold",color=CRIMSON)
ax.text(RX,2.08,"≈ 0.8 for this CpG",ha="left",va="center",fontsize=9.5,color=TEXT)
ax.text(RX,1.65,"+100 keeps β stable when\nboth signals are weak",ha="left",va="center",fontsize=8.8,color=TEXT2,
        linespacing=1.4)

# Infinium I note
ax.add_patch(mp.FancyBboxPatch((0.0,-0.95),11.95,1.0,boxstyle="round,pad=0.02,rounding_size=0.08",
                               fc=FAINT,ec="none"))
ax.text(0.15,-0.15,"Infinium I (a minority of probes)",ha="left",va="center",fontsize=9.5,fontweight="semibold",color=TEXT)
ax.text(0.15,-0.58,"two bead types per CpG — an M probe and a U probe whose 3′ base matches C or T.\n"
        "Only the matching probe is extended; both are read in the same colour channel.",ha="left",va="center",linespacing=1.5,
        fontsize=8.8,color=TEXT2)
ax.set_xlim(-0.1,12.1); ax.set_ylim(-1.05,5.5); ax.axis("off")
header(ax,"DNA methylation · Illumina Infinium array (EPIC)",
       "The array reads each CpG as two colours of one extra base",dy=6)
fig.savefig(OUT+"array_principle.png"); print("ok g20 array")
