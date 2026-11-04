from gallery import *

fig,(ax,bx,cx)=plt.subplots(3,1,figsize=(11.5,10.4),
                            gridspec_kw={"height_ratios":[0.95,1.2,1.0],"hspace":0.3})

# L1 parts in kb (L1RP-like consensus): 5'UTR, ORF1, ORF2 (EN, RT domains), 3'UTR, poly-A
PARTS=[(0.00,0.91,tint(0.30),"5′UTR"),(0.91,1.93,tint(0.50),"ORF1"),
       (1.99,5.82,tint(0.50),"ORF2"),(5.82,6.03,tint(0.30),""),(6.03,6.22,tint(0.30),"")]
DOMS=[(2.00,2.72,"EN"),(3.48,4.32,"RT")]
H=0.5

def l1(axx,y,x0=0.0,s=1.0,keep=(0,9),faded=False,labels=True):
    """Draw an L1 element at height y; keep=(a,b) in kb shows only that stretch."""
    a,b=keep
    for p,q,c,name in PARTS:
        p2,q2=max(p,a),min(q,b)
        if p2>=q2: continue
        axx.add_patch(mp.Rectangle((x0+p2*s,y),(q2-p2)*s,H,fc=FAINT if faded else c,ec="white",lw=1.2,zorder=2))
        if labels and name and p>=a:
            axx.text(x0+(p+q)/2*s,y+H+0.12,name,ha="center",va="bottom",fontsize=11,
                     fontweight="semibold",color=TEXT)
    for p,q,name in DOMS:
        if p<a or faded: continue
        axx.add_patch(mp.Rectangle((x0+p*s,y+0.07),(q-p)*s,H-0.14,fc=CRIMSON,ec="none",zorder=3))
        axx.text(x0+(p+q)/2*s,y+H/2,name,ha="center",va="center",fontsize=10.5,
                 fontweight="bold",color="white",zorder=4)

# ① anatomy ------------------------------------------------------------------
ax.plot([-0.55,0],[0.25,0.25],color=LINE,lw=2); ax.plot([6.22,6.8],[0.25,0.25],color=LINE,lw=2)
ax.text(-0.6,0.25,"genome",ha="right",va="center",fontsize=9.5,color=TEXT2)
l1(ax,0)
ax.text(6.02,H+0.12,"3′UTR · A(n)",ha="center",va="bottom",fontsize=9.5,color=TEXT2)
roles=[(0.45,0,"internal promoter","makes an RNA copy of itself"),
       (1.42,1,"ORF1p","RNA-binding chaperone"),
       (2.36,0,"endonuclease (EN)","cuts the target DNA"),
       (3.90,1,"reverse transcriptase (RT)","copies RNA back into DNA"),
       (6.13,0,"poly-A tail","where copying starts")]
for x,lv,t1,t2 in roles:
    y=-0.36-lv*0.72
    ax.plot([x,x],[-0.05,y+0.08],color=LINE,lw=0.9)
    ax.text(x,y,t1,ha="center",va="top",fontsize=10,fontweight="semibold",color=TEXT)
    ax.text(x,y-0.3,t2,ha="center",va="top",fontsize=9.5,color=TEXT2)
ax.plot([-1.2,-0.2],[-1.25,-1.25],color=TEXT2,lw=1.2)
ax.text(-0.7,-1.33,"1 kb",ha="center",va="top",fontsize=9,color=TEXT2)
ax.set_xlim(-1.3,7.2); ax.set_ylim(-1.75,1.0); ax.axis("off")
header(ax,"Transposons · how an L1 copies itself",
       "A transposon carries the blueprint for its own copy machine",x=0.0,dy=26)
ax.set_title("① Anatomy of an active L1 (≈ 6 kb)",loc="left",fontsize=12.5,
             fontweight="semibold",color=TEXT,pad=8)

# ② copy & paste ---------------------------------------------------------------
CW=2.75                                   # column width
def dna(x0,x1,y,nick=None):
    for yy in (y,y-0.16):
        if nick is not None and yy<y:
            bx.plot([x0,nick-0.04],[yy,yy],color=TEXT2,lw=1.6)
            bx.plot([nick+0.04,x1],[yy,yy],color=TEXT2,lw=1.6)
        else:
            bx.plot([x0,x1],[yy,yy],color=TEXT2,lw=1.6)
def rna(x0,x1,y):
    bx.plot([x0,x1],[y,y],color=CRIMSON,lw=1.8,solid_capstyle="round")
    bx.text(x1+0.05,y,"AAA",ha="left",va="center",fontsize=8.5,color=CRIMSON)
def blob(x,y,w,h,c,label=None,lc=TEXT2,dy=0.22):
    bx.add_patch(mp.Ellipse((x,y),w,h,fc=c,ec="none",zorder=3))
    if label: bx.text(x,y+dy,label,ha="center",va="bottom",fontsize=9,color=lc)

steps=[("Transcribe","L1 is copied into RNA —","the original stays put"),
       ("Make proteins","ORF1p and ORF2p are made","and bind their own RNA"),
       ("Nick & copy","ORF2p cuts a new site and","copies the RNA into DNA there"),
       ("New copy","a new L1 sits in a new place,","flanked by short duplications")]
for i,(t,c1,c2) in enumerate(steps):
    x0=i*CW
    bx.add_patch(mp.Circle((x0+0.2,2.55),0.15,fc=CRIMSON,ec="none"))
    bx.text(x0+0.2,2.55,str(i+1),ha="center",va="center",fontsize=10,fontweight="bold",color="white")
    bx.text(x0+0.45,2.55,t,ha="left",va="center",fontsize=11.5,fontweight="semibold",color=TEXT)
    bx.text(x0+0.05,0.12,c1,ha="left",va="top",fontsize=9.5,color=TEXT2)
    bx.text(x0+0.05,-0.18,c2,ha="left",va="top",fontsize=9.5,color=TEXT2)
    if i: bx.text(x0-0.12,1.35,"›",ha="center",va="center",fontsize=26,color=LINE)
    bx.text(x0+0.05,2.18,"nucleus" if i!=1 else "cytoplasm",ha="left",va="center",
            fontsize=8.5,color=TEXT2,style="italic")

# 1 transcribe
x0=0
dna(x0+0.05,x0+2.5,1.75)
bx.add_patch(mp.Rectangle((x0+0.7,1.55),1.15,0.24,fc=tint(0.5),ec="none",zorder=2))
bx.text(x0+1.27,1.92,"L1",ha="center",va="bottom",fontsize=9.5,fontweight="semibold",color=TEXT)
bx.annotate("",(x0+1.27,0.98),(x0+1.27,1.45),arrowprops=dict(arrowstyle="-|>",color=TEXT2,lw=1.2))
rna(x0+0.7,x0+1.85,0.82); bx.text(x0+0.6,0.82,"RNA",ha="right",va="center",fontsize=9,color=CRIMSON)
# 2 proteins on RNA
x0=CW
rna(x0+0.35,x0+2.0,1.05)
for k in range(4): blob(x0+0.55+k*0.3,1.17,0.24,0.2,tint(0.45))
bx.text(x0+1.0,1.38,"ORF1p",ha="center",va="bottom",fontsize=9,color=TEXT2)
blob(x0+1.82,1.24,0.42,0.34,CRIMSON,"ORF2p",dy=0.2)
# 3 nick & reverse transcribe (target-primed)
x0=2*CW; nk=x0+1.95
dna(x0+0.05,x0+2.5,1.75,nick=nk)
blob(nk,1.67,0.34,0.36,CRIMSON); bx.text(nk,1.97,"EN cuts",ha="center",va="bottom",fontsize=9,color=CRIMSON)
rna(x0+0.45,nk-0.18,0.95); bx.text(x0+0.35,0.95,"RNA",ha="right",va="center",fontsize=9,color=CRIMSON)
bx.plot([nk,nk,nk-0.1],[1.5,1.22,1.22],color=CRIMSON,lw=3.2,solid_joinstyle="round")
bx.annotate("",(x0+0.7,1.22),(nk-0.05,1.22),
            arrowprops=dict(arrowstyle="-|>",color=CRIMSON,lw=3.2,mutation_scale=14))
bx.text(x0+1.2,0.78,"RT copies RNA → DNA",ha="center",va="top",fontsize=9,color=CRIMSON,fontweight="semibold")
# 4 new copy with target-site duplications
x0=3*CW
dna(x0+0.05,x0+2.5,1.75)
bx.add_patch(mp.Rectangle((x0+0.85,1.55),0.95,0.24,fc=tint(0.5),ec="none",zorder=2))
for xx in (x0+0.62,x0+1.8):
    bx.add_patch(mp.Rectangle((xx,1.55),0.23,0.24,fc="white",ec=CRIMSON,lw=1,zorder=2))
bx.text(x0+1.32,1.92,"new L1",ha="center",va="bottom",fontsize=9.5,fontweight="semibold",color=TEXT)
bx.text(x0+0.73,1.42,"TSD",ha="center",va="top",fontsize=8.5,color=TEXT2)
bx.text(x0+1.91,1.42,"TSD",ha="center",va="top",fontsize=8.5,color=TEXT2)
bx.plot([0,4*CW-0.2],[-0.68,-0.68],color="#E3E6EB",lw=0.9)
bx.text(0.05,-0.82,"Copy & paste: the original never leaves, so copies pile up (L1 ≈ 20% of our genome).  "
        "Alu (≈ 300 bp) has no genes of its own — it borrows L1's ORF2p.",
        ha="left",va="top",fontsize=9.5,color=TEXT)
bx.set_xlim(-0.05,4*CW); bx.set_ylim(-1.05,2.75); bx.axis("off")
bx.set_title("② Copy & paste through an RNA intermediate (retrotransposition)",loc="left",
             fontsize=12.5,fontweight="semibold",color=TEXT,pad=8)

# ③ fossils ----------------------------------------------------------------------
S=1.0; X0=0.0
rows=[("intact","≈ 80–100 per genome",dict(),
       "can still jump",True),
      ("5′-truncated","most insertions",dict(keep=(5.0,6.22),labels=False),
       "copying starts at the poly-A end, often stops early → front lost: no promoter, no ORFs",False),
      ("full-length, old","mutated over time",dict(faded=True,labels=False),
       "nothing keeps it intact → stop codons, frameshifts pile up → no working ORF2p",False)]
for i,(lab,sub,kw,note,good) in enumerate(rows):
    y=-i*1.05
    cx.text(-0.25,y+H/2+0.1,lab,ha="right",va="center",fontsize=11,fontweight="semibold",color=TEXT)
    cx.text(-0.25,y+H/2-0.17,sub,ha="right",va="center",fontsize=9.5,color=TEXT2)
    if "keep" in kw:
        cx.add_patch(mp.Rectangle((X0,y),5.0,H,fc="white",ec=LINE,lw=1,ls=(0,(3,2)),zorder=1))
        cx.text(2.5,y+H/2,"never copied",ha="center",va="center",fontsize=9.5,color=TEXT2)
    l1(cx,y,X0,S,**kw)
    if kw.get("faded"):
        for xm,mk in [(0.55,"×"),(1.35,"×"),(2.6,"×"),(3.2,"stop"),(4.05,"×"),(4.7,"fs"),(5.35,"×")]:
            cx.text(xm,y+H/2,mk,ha="center",va="center",fontsize=12 if mk=="×" else 9,
                    fontweight="bold",color=CRIMSON)
    if good:
        cx.text(6.55,y+H/2,note,ha="left",va="center",fontsize=10.5,color=CRIMSON,fontweight="bold")
    else:
        w=note.split(" → ")
        cx.text(6.55,y+H/2+0.14,w[0]+" →",ha="left",va="center",fontsize=9.5,color=TEXT2)
        cx.text(6.55,y+H/2-0.16," → ".join(w[1:]),ha="left",va="center",fontsize=9.5,color=TEXT,
                fontweight="semibold")
cx.text(-2.2,-2.75,"≈ 500,000 L1 copies in the human genome, but only ≈ 80–100 can still jump — "
        "the rest are fossils (they can also be silenced by DNA methylation).",
        ha="left",va="center",fontsize=9.5,color=TEXT)
cx.set_xlim(-2.2,11.6); cx.set_ylim(-2.95,0.95); cx.axis("off")
cx.set_title("③ Why most copies are fossils",loc="left",fontsize=12.5,fontweight="semibold",
             color=TEXT,pad=8,x=-0.0)
fig.savefig(OUT+"transposon.png"); print("ok f1 transposon")
