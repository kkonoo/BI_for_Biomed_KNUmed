from gallery import *

fig,ax=plt.subplots(figsize=(11.5,7.2))
CW,CH=0.42,0.34                         # cell width / height
MONO="DejaVu Sans Mono"
REF=list("ATCGACTCGACA")                # CpG1 = 2–3 (methylated), CpG2 = 7–8 (unmethylated); 5, 10 = non-CpG C
M1,M2,NC=2,7,(5,10)

def cell(x,y,ch,kind="plain",w=CW,h=CH,fs=12):
    fc,tc,ec={"plain":(FAINT,TEXT,"none"),"meth":(CRIMSON,"white","none"),
              "conv":(tint(0.18),CRIMSON,"none"),"t":(FAINT,"#9AA1AB","none"),
              "tconv":("white",TEXT2,LINE)}[kind]
    ax.add_patch(mp.Rectangle((x,y),w-0.04,h,fc=fc,ec=ec,lw=0.8))
    ax.text(x+(w-0.04)/2,y+h/2,ch,ha="center",va="center",fontsize=fs,family=MONO,color=tc,
            fontweight="bold" if kind in ("meth","conv") else "normal")

def row(x0,y,seq,kinds):
    for i,(ch,k) in enumerate(zip(seq,kinds)): cell(x0+i*CW,y,ch,k)

X0=1.75
steps=[("Genomic DNA","5mC at CpG 1, plain C at CpG 2",REF,
        ["meth" if i==M1 else "plain" for i in range(12)]),
       ("① Bisulfite","unmethylated C → U;  5mC stays C",
        [("U" if (c=="C" and i!=M1) else c) for i,c in enumerate(REF)],
        ["meth" if i==M1 else ("conv" if (REF[i]=="C") else "plain") for i in range(12)]),
       ("② PCR","U is copied as T",
        [("T" if (c=="C" and i!=M1) else c) for i,c in enumerate(REF)],
        ["meth" if i==M1 else ("conv" if (REF[i]=="C") else "plain") for i in range(12)])]
for k,(t,sub,seq,kinds) in enumerate(steps):
    y=2.95-k*1.02
    ax.text(X0-0.2,y+CH/2+0.08,t,ha="right",va="center",fontsize=11,fontweight="semibold",color=TEXT)
    ax.text(X0-0.2,y+CH/2-0.17,sub,ha="right",va="center",fontsize=8.8,color=TEXT2)
    row(X0,y,seq,kinds)
    if k: ax.annotate("",(X0+12*CW+0.2,y+CH/2+0.05),(X0+12*CW+0.2,y+1.02+CH/2-0.05),
                      arrowprops=dict(arrowstyle="-|>",color=LINE,lw=1.2,mutation_scale=11))
# CpG labels on the genomic row
yg=2.95
for i,lab,filled in [(M1,"CpG 1",True),(M2,"CpG 2",False)]:
    xc=X0+i*CW+(CW-0.04)/2+0.02
    ax.plot([xc-0.02,xc-0.02],[yg+CH,yg+CH+0.16],color=TEXT2,lw=1)
    ax.add_patch(mp.Circle((xc-0.02,yg+CH+0.24),0.075,fc=CRIMSON if filled else "white",ec=CRIMSON,lw=1.3))
    ax.text(xc-0.02,yg+CH+0.4,lab,ha="center",va="bottom",fontsize=9.5,fontweight="semibold",color=TEXT)
for i in NC:
    ax.text(X0+i*CW+CW/2,yg-0.05,"non-CpG",ha="center",va="top",fontsize=7.5,color=TEXT2)

# ③ reads aligned to the reference
RX=7.75; ry=3.0
ax.text(RX,ry+0.62,"③ Align reads to the reference",ha="left",va="center",fontsize=11,fontweight="semibold",color=TEXT)
ax.text(RX,ry+0.36,"compare each C position: C = was methylated, T = was not",ha="left",va="center",fontsize=8.8,color=TEXT2)
ax.text(RX-0.12,ry-0.12+CH/2,"ref",ha="right",va="center",fontsize=9,color=TEXT2)
for i,c in enumerate(REF): cell(RX+i*0.3,ry-0.12,c,"plain",w=0.3,h=0.28,fs=9.5)
reads=10; rh=0.2
c1=[1]*9+[0]; c2=[0]*4+[1]+[0]*5
for r in range(reads):
    y=ry-0.5-r*rh
    a,b=0,12
    ax.add_patch(mp.Rectangle((RX+a*0.3,y+0.03),(b-a)*0.3-0.04,rh-0.06,fc=FAINT,ec="none"))
    for i in range(a,b):
        if REF[i]!="C": continue
        if i==M1: ch,k=("C","meth") if c1[r] else ("T","tconv")
        elif i==M2: ch,k=("C","meth") if c2[r] else ("T","tconv")
        else: ch,k="T","t"
        fc,tc={"meth":(CRIMSON,"white"),"tconv":("white",TEXT2),"t":(FAINT,"#9AA1AB")}[k]
        ax.add_patch(mp.Rectangle((RX+i*0.3,y+0.01),0.26,rh-0.02,fc=fc,ec=LINE if k=="tconv" else "none",lw=0.6))
        ax.text(RX+i*0.3+0.13,y+rh/2,ch,ha="center",va="center",fontsize=7.5,family=MONO,color=tc,
                fontweight="bold" if k=="meth" else "normal")
ax.text(RX+12*0.3+0.15,ry-0.5-4.5*rh+rh/2,"10 reads",ha="left",va="center",fontsize=9,color=TEXT2)

# ④ per-CpG counts → β
by=ry-0.62-reads*rh-0.25
ax.text(RX,by,"④ Count per CpG",ha="left",va="center",fontsize=11,fontweight="semibold",color=TEXT)
for j,(lab,cc,tt,beta) in enumerate([("CpG 1",9,1,0.9),("CpG 2",1,9,0.1)]):
    y=by-0.36-j*0.34
    ax.text(RX,y,f"{lab}:  {cc} C / {cc+tt} reads",ha="left",va="center",fontsize=10,color=TEXT)
    ax.text(RX+2.25,y,f"β = {beta:.1f}",ha="left",va="center",fontsize=11,fontweight="bold",
            color=CRIMSON if beta>0.5 else TEXT2)
ax.text(RX,by-1.12,"β = C / (C + T)",ha="left",va="center",fontsize=10,fontweight="semibold",color=TEXT)
ax.text(RX+1.55,by-1.12,"  non-CpG C always reads T → checks that conversion worked",ha="left",va="center",
        fontsize=8.8,color=TEXT2)

# legend + caveat
for k,(kind,ch,t) in enumerate([("meth","C","5mC: protected, still reads C"),("conv","U","unmethylated C: converted to U, then T"),
                              ("plain","A","other bases: unchanged")]):
    ly=0.15-k*0.38
    cell(X0,ly,ch,kind,w=0.34,h=0.28,fs=9.5)
    ax.text(X0+0.42,ly+0.14,t,ha="left",va="center",fontsize=9,color=TEXT2)
ax.text(X0-2.6,-1.55,"Caveat: 5-hydroxymethylcytosine (5hmC) is also protected and reads as C — "
        "bisulfite cannot tell 5mC from 5hmC.",ha="left",va="center",fontsize=8.8,color=TEXT2)
ax.set_xlim(-1.0,12.6); ax.set_ylim(-1.75,4.0); ax.axis("off")
header(ax,"DNA methylation · bisulfite sequencing","Bisulfite turns a chemical mark into a sequence difference",
       x=0.0,dy=-8)
fig.savefig(OUT+"bisulfite_principle.png"); print("ok g19 bisulfite")
