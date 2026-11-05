from gallery import *
from matplotlib.colors import LinearSegmentedColormap
rng=np.random.default_rng(17)
CM=LinearSegmentedColormap.from_list("hic",["#FFFFFF",tint(0.25),tint(0.6),CRIMSON,"#5E1220"])
DARK="#5A5F66"

FW,FH=11.5,6.6
fig=plt.figure(figsize=(FW,FH))
W=0.29; X0=[0.0,0.355,0.71]; MH=W*FW/FH; MY=0.355
def axes(i):
    m=fig.add_axes([X0[i],MY,W,MH]); t=fig.add_axes([X0[i],MY-0.085,W,0.065])
    c=fig.add_axes([X0[i],0.0,W,0.19]); return m,t,c
def decay(n,a=1.0):
    i=np.arange(n); d=np.abs(i[:,None]-i[None,:]); return (d+1.0)**(-a)
def show(m,M,title,scale):
    M=rng.poisson(M*400)/400+1e-4
    L=np.log10(M)
    m.imshow(L,cmap=CM,interpolation="nearest",vmin=np.percentile(L,30),vmax=np.percentile(L,99.7))
    m.set_xticks([]); m.set_yticks([])
    for s in m.spines.values(): s.set_color(LINE); s.set_visible(True)
    m.set_title(title,loc="left",fontsize=12.5,fontweight="semibold",color=TEXT,pad=22)
    m.text(0,1.025,scale,transform=m.transAxes,ha="left",va="bottom",fontsize=9.5,color=TEXT2)

# (a) A/B compartments -----------------------------------------------------------------
n=120; seg=[0,14,30,41,58,70,88,101,120]
c=np.zeros(n)
for k in range(len(seg)-1): c[seg[k]:seg[k+1]]=1 if k%2==0 else -1
ev=np.convolve(c,np.ones(5)/5,mode="same")*(0.75+0.25*rng.random(n))
M=decay(n,0.75)*(1+0.85*np.outer(c,c)*0.5+0.0)
m,t,cx=axes(0); show(m,M,"(a) A/B compartments","multi-Mb scale · plaid pattern")
t.bar(np.arange(n),ev,width=1.0,color=[CRIMSON if v>0 else "#9AA1AB" for v in ev],lw=0)
t.axhline(0,color=LINE,lw=0.8); t.set_xlim(-0.5,n-0.5); t.set_ylim(-1.15,1.15); t.axis("off")
t.text(-2,0.55,"A",ha="right",va="center",fontsize=10,fontweight="bold",color=CRIMSON)
t.text(-2,-0.55,"B",ha="right",va="center",fontsize=10,fontweight="bold",color=DARK)
t.text(n/2,-1.35,"first eigenvector (PC1)",ha="center",va="top",fontsize=8.8,color=TEXT2)
cx.text(0.5,0.7,"A = active, gene-rich, open\nB = inactive, gene-poor, closed",ha="center",va="center",
        fontsize=9.5,color=TEXT,transform=cx.transAxes,linespacing=1.5)
cx.text(0.5,0.22,"A contacts A, B contacts B — even far apart",ha="center",va="center",fontsize=9.5,
        color=CRIMSON,fontweight="semibold",transform=cx.transAxes)
cx.axis("off")
header(m,"3D genome · reading a Hi-C contact matrix","Three features at three scales",dy=40,x1=1/W)

# (b) TADs ----------------------------------------------------------------------------------
n=100; bd=[0,22,47,71,100]
lab=np.zeros(n,int)
for k in range(len(bd)-1): lab[bd[k]:bd[k+1]]=k
same=(lab[:,None]==lab[None,:])
M=decay(n,0.9)*(1+2.2*same)
m,t,cx=axes(1); show(m,M,"(b) TADs","~100 kb – 1 Mb · blocks on the diagonal")
for b in bd[1:-1]:
    m.plot([b-0.5,b-0.5],[b-6,b+6],color=TEXT,lw=0,alpha=0)
ins=np.array([M[i-6:i,i+1:i+7].mean() if 6<=i<n-7 else np.nan for i in range(n)])
ins=np.log2(ins/np.nanmean(ins))
t.plot(np.arange(n),ins,color=CRIMSON,lw=1.4); t.set_xlim(-0.5,n-0.5); t.axis("off")
for b in bd[1:-1]: t.axvline(b-0.5,color=LINE,lw=0.8,ls=(0,(2,2)))
t.text(n/2,np.nanmin(ins)-0.35,"insulation score: dips = boundaries",ha="center",va="top",fontsize=8.8,color=TEXT2)
t.set_ylim(np.nanmin(ins)-0.3,np.nanmax(ins)+0.2)
# cartoon: domains as loops of chromatin, CTCF at the boundaries
cx.set_xlim(0,100); cx.set_ylim(0,1); cx.axis("off")
cx.plot([0,100],[0.32,0.32],color=TEXT2,lw=1.2)
for (a,b) in zip(bd[:-1],bd[1:]):
    w=(b-a)*0.82; xc=(a+b)/2
    cx.add_patch(mp.Ellipse((xc,0.58),w,0.5,fc=tint(0.12),ec=tint(0.6),lw=1.1))
for b in bd[1:-1]:
    cx.add_patch(mp.Polygon([(b-1.8,0.2),(b+1.8,0.2),(b,0.36)],closed=True,fc=CRIMSON,ec="none",zorder=3))
cx.text(bd[1],0.03,"CTCF",ha="center",va="bottom",fontsize=9,color=CRIMSON,fontweight="semibold")
cx.text(bd[2],0.03,"CTCF",ha="center",va="bottom",fontsize=9,color=CRIMSON,fontweight="semibold")
cx.text(bd[3],0.03,"CTCF",ha="center",va="bottom",fontsize=9,color=CRIMSON,fontweight="semibold")
cx.text(50,0.98,"contacts stay inside a domain",ha="center",va="top",fontsize=9.5,color=TEXT)

# (c) loops -------------------------------------------------------------------------------------
n=100; a0,b0=22,74
lab=np.zeros(n,int); lab[a0:b0+1]=1
same=(lab[:,None]==lab[None,:])&(lab[:,None]==1)
M=decay(n,0.95)*(1+1.2*same)
i=np.arange(n)
dot=6.0*np.exp(-(((i[:,None]-a0)**2+(i[None,:]-b0)**2)/(2*1.3**2)))
M=M+(dot+dot.T)*M.max()*0.02
m,t,cx=axes(2); show(m,M,"(c) Loops","~10–100 kb · a dot off the diagonal")
m.add_patch(mp.Circle((b0,a0),4.2,fill=False,ec=TEXT,lw=1.3))
m.annotate("loop",(b0+2.8,a0-3.5),(b0+12,a0-14),fontsize=10,color=TEXT,fontweight="semibold",
           ha="left",va="center",arrowprops=dict(arrowstyle="-",color=TEXT,lw=1))
t.set_xlim(-0.5,n-0.5); t.set_ylim(0,1); t.axis("off")
t.plot([0,n],[0.5,0.5],color=LINE,lw=0.9)
for xa,d in [(a0,1),(b0,-1)]:
    t.add_patch(mp.FancyArrow(xa-d*2.5,0.5,d*4.5,0,width=0.32,head_width=0.6,head_length=2.2,
                              length_includes_head=True,fc=CRIMSON,ec="none"))
t.text(n/2,-0.25,"loop anchors: two CTCF sites, facing each other",ha="center",va="top",fontsize=8.8,color=TEXT2)
# cartoon: a chromatin loop held by a cohesin ring
cx.set_xlim(0,100); cx.set_ylim(0,1); cx.axis("off")
th=np.radians(np.linspace(-68,248,200))
cx.plot(50+13*np.cos(th),0.63+0.29*np.sin(th),color=TEXT2,lw=1.3)
lx,ly=50+13*np.cos(th[0]),0.63+0.29*np.sin(th[0])
cx.plot([8,100-lx],[0.3,ly],color=TEXT2,lw=1.3); cx.plot([lx,92],[ly,0.3],color=TEXT2,lw=1.3)
cx.add_patch(mp.Rectangle((36,0.29),6.5,0.1,fc=tint(0.5),ec="none",zorder=3,angle=0))
cx.add_patch(mp.Rectangle((57.5,0.29),6.5,0.1,fc=TEXT,ec="none",zorder=3))
cx.add_patch(mp.Ellipse((50,0.36),11,0.26,fill=False,ec=CRIMSON,lw=2.6,zorder=4))
cx.text(50,0.66,"cohesin ring\nholds the base",ha="center",va="center",fontsize=8.8,color=CRIMSON)
cx.text(37,0.16,"enhancer",ha="center",va="center",fontsize=9.5,color=TEXT)
cx.text(63,0.16,"promoter",ha="center",va="center",fontsize=9.5,color=TEXT)
cx.text(50,-0.02,"two distant sites touch; the DNA between loops out",ha="center",va="bottom",fontsize=9,color=TEXT2)

fig.text(0.0,-0.05,"Simulated matrices (schematic). Each pixel = how often two genomic bins were ligated together; darker = more contacts.\n"
         "TAD boundaries are largely shared between cell types; loops are more cell-type specific.",
         ha="left",va="top",fontsize=8.8,color=TEXT2,linespacing=1.5)
fig.savefig(OUT+"hic_features.png"); print("ok g17 hic features")
