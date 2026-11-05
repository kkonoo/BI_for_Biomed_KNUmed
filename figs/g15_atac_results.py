from gallery import *
rng=np.random.default_rng(15)

fig=plt.figure(figsize=(11.5,5.6))
ax=fig.add_axes([0.0,0.36,0.56,0.47])      # (a) ATAC track
gx=fig.add_axes([0.0,0.06,0.56,0.27])      # (a) gene model
cx=fig.add_axes([0.65,0.52,0.35,0.31])     # (b) cartoon
fx=fig.add_axes([0.65,0.08,0.35,0.36])     # (b) aggregate footprint

# (a) open chromatin = peaks --------------------------------------------------------
x=np.linspace(0,60,3000)
els=[(9.0,0.55,0.62,"enhancer"),(22.0,0.35,1.0,"promoter (TSS)"),(37.0,0.5,0.7,"enhancer")]
y=0.03+0.025*rng.random(x.size)
for c,s,h,_ in els: y+=h*np.exp(-0.5*((x-c)/s)**2)
for c,s,h,lab in els:
    ax.axvspan(c-1.6,c+1.6,color=FAINT,zorder=0); gx.axvspan(c-1.6,c+1.6,color=FAINT,zorder=0)
ax.fill_between(x,0,y,color=tint(0.55),lw=0,zorder=2)
ax.plot(x,y,color=CRIMSON,lw=0.8,zorder=3)
ax.text(0.3,1.22,"ATAC-seq coverage",ha="left",va="top",fontsize=10.5,fontweight="semibold",color=TEXT)
for c,s,h,lab in els:
    ax.text(c,h+0.13,"peak",ha="center",va="bottom",fontsize=9,color=CRIMSON,fontweight="semibold")
ax.text(50.5,0.22,"closed chromatin:\nno insertions, no reads",ha="center",va="center",fontsize=9,color=TEXT2)
ax.set_xlim(0,60); ax.set_ylim(0,1.25); ax.axis("off")
ax.set_title("(a) Open chromatin = peaks",loc="left",fontsize=12.5,fontweight="semibold",color=TEXT,pad=6)
header(ax,"Chromatin accessibility · what ATAC-seq gives you",
       "Peaks mark open elements; footprints mark bound factors",dy=22,x1=1/0.56)

# gene model
gx.plot([22,48],[0.62,0.62],color=TEXT2,lw=1.2)
for a,b in [(22,23.2),(28,28.9),(33,33.8),(41,41.9),(46.5,48)]:
    gx.add_patch(mp.Rectangle((a,0.5),b-a,0.24,fc=TEXT,ec="none"))
gx.annotate("",(24.6,0.98),(22,0.98),arrowprops=dict(arrowstyle="-|>",color=TEXT,lw=1.1,mutation_scale=10))
gx.plot([22,22],[0.74,0.98],color=TEXT,lw=1.1)
gx.text(35,0.92,"GENE",ha="center",va="bottom",fontsize=10,style="italic",color=TEXT)
for c,s,h,lab in els:
    gx.text(c,0.22,lab,ha="center",va="center",fontsize=9.5,fontweight="semibold",
            color=CRIMSON if "promoter" in lab else TEXT)
gx.plot([0,60],[0.03,0.03],color=LINE,lw=0.9)
gx.plot([49,59],[-0.12,-0.12],color=TEXT2,lw=1.1)
gx.text(54,-0.2,"10 kb",ha="center",va="top",fontsize=9,color=TEXT2)
gx.set_xlim(0,60); gx.set_ylim(-0.4,1.2); gx.axis("off")

# (b) cartoon: TF bound on its motif, Tn5 cuts only beside it ----------------------
for yy in (0.40,0.30): cx.plot([0,10],[yy,yy],color=TEXT2,lw=1.6)
cx.add_patch(mp.Rectangle((4.3,0.28),1.4,0.14,fc=tint(0.35),ec="none",zorder=2))
cx.text(5,0.18,"motif",ha="center",va="top",fontsize=9,color=TEXT2)
cx.add_patch(mp.Ellipse((5,0.62),1.9,0.42,fc=CRIMSON,ec="none",zorder=3))
cx.text(5,0.62,"TF",ha="center",va="center",fontsize=10.5,fontweight="bold",color="white",zorder=4)
for xc in (1.6,2.6,3.5,6.5,7.4,8.4):
    cx.annotate("",(xc,0.44),(xc,0.80),arrowprops=dict(arrowstyle="-|>",color=TEXT2,lw=1.1,mutation_scale=9))
cx.text(2.55,0.86,"Tn5 cuts here",ha="center",va="bottom",fontsize=9,color=TEXT2)
cx.text(7.45,0.86,"and here",ha="center",va="bottom",fontsize=9,color=TEXT2)
cx.text(5,0.98,"protected",ha="center",va="bottom",fontsize=9,color=CRIMSON,fontweight="semibold")
cx.set_xlim(0,10); cx.set_ylim(0.05,1.12); cx.axis("off")
cx.set_title("(b) Footprint: a bound factor protects its own site",loc="left",fontsize=12.5,
             fontweight="semibold",color=TEXT,pad=6)

# aggregate Tn5 insertions around many motif instances
p=np.arange(-100,101)
f=1+0.75*np.exp(-0.5*((np.abs(p)-24)/11)**2)-0.55*np.exp(-0.5*(p/6.5)**2)
f=f+0.04*rng.standard_normal(p.size)
fx.axvspan(-9,9,color=FAINT,zorder=0)
fx.plot(p,f,color=CRIMSON,lw=1.8)
fx.axhline(1,color=LINE,lw=0.9,ls=(0,(3,2)))
fx.annotate("footprint",(0,0.55),(42,0.55),ha="left",va="center",fontsize=9.5,color=CRIMSON,
            fontweight="semibold",arrowprops=dict(arrowstyle="-|>",color=CRIMSON,lw=1,mutation_scale=9))
fx.set_xlim(-100,100); fx.set_ylim(0.3,2.05)
fx.set_xlabel("distance from motif centre (bp)",color=TEXT2,fontsize=9.5)
fx.set_ylabel("Tn5 insertions\n(relative)",color=TEXT2,fontsize=9.5)
fx.tick_params(colors=TEXT2,labelsize=9); fx.set_yticks([0.5,1,1.5,2])
for s in ("left","bottom"): fx.spines[s].set_color(LINE)
fx.text(-97,1.98,"average over thousands of motif sites",ha="left",va="top",fontsize=8.8,color=TEXT2)

fig.text(0.65,-0.06,"Caveat: Tn5's own sequence preference can create footprint-like dips, and factors\n"
         "that bind only briefly leave no footprint at all.  Schematic.",ha="left",va="top",
         fontsize=8.8,color=TEXT2)
fig.savefig(OUT+"atac_results.png"); print("ok g15 atac results")
