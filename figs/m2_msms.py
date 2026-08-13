from style import *
import matplotlib.patches as mp
rng=np.random.default_rng(11)
fig=plt.figure(figsize=(12.6,6.4))
gs=fig.add_gridspec(2,1,height_ratios=[0.85,1.5],hspace=0.42)

# --- workflow strip ---
ax=fig.add_subplot(gs[0]); ax.axis("off")
steps=[("Protein\nextract","#9AA5B1"),("Trypsin\ndigest",TEAL),("LC\nseparation",BLUE),
       ("MS1\nsurvey","#C9A227"),("Fragment\n(HCD)","#E08A3C"),("MS2\nspectrum",PURPLE),
       ("Database\nsearch",RED),("Protein\nlist",NAVY)]
w=0.108; gap=0.0175
for i,(t,c) in enumerate(steps):
    x=i*(w+gap)
    ax.add_patch(mp.FancyBboxPatch((x,0.42),w,0.44,boxstyle="round,pad=0.006,rounding_size=0.02",
                 fc=c,ec="none",alpha=0.85))
    ax.text(x+w/2,0.64,t,ha="center",va="center",color="white",fontsize=9.6,
            fontweight="bold",linespacing=1.3)
    if i<len(steps)-1:
        ax.annotate("",xy=(x+w+gap*0.92,0.64),xytext=(x+w+gap*0.08,0.64),
                    arrowprops=dict(arrowstyle="-|>",color="#999",lw=1.4))
ax.text(0.0,0.16,"⚠️  You never observe a protein. You observe peptides, and infer the protein (§3.2).",
        fontsize=10.2,color="#8A3A33")
ax.set_xlim(-0.01,1.0); ax.set_ylim(0,1)
ax.set_title("The bottom-up proteomics workflow",loc="left",x=0,y=1.0,fontsize=13)

# --- MS2 spectrum ---
ax=fig.add_subplot(gs[1])
pep="P E P T I D E K".split()
bion=[98,227,324,425,538,653,768]
yion=[147,262,377,490,591,720,849]
for mz,h,c in [(m,rng.uniform(18,95),BLUE) for m in bion]:
    ax.vlines(mz,0,h,color=c,lw=2.0)
    ax.text(mz,h+2.5,f"b{bion.index(mz)+1}",ha="center",fontsize=8.4,color=c)
for mz,h,c in [(m,rng.uniform(25,100),RED) for m in yion]:
    ax.vlines(mz,0,h,color=c,lw=2.0)
    ax.text(mz,h+2.5,f"y{yion.index(mz)+1}",ha="center",fontsize=8.4,color=c)
noise_mz=rng.uniform(80,950,240); noise_h=rng.gamma(1.1,2.4,240)
ax.vlines(noise_mz,0,noise_h,color="#CFCFCF",lw=0.9,zorder=0)
ax.set_xlabel("m/z"); ax.set_ylabel("relative intensity")
ax.set_xlim(60,980); ax.set_ylim(0,118)
ax.set_title("An MS2 spectrum: the peptide's sequence, read from the gaps",loc="left",fontsize=12.5)
ax.text(0.985,0.95,"blue = b ions (N-terminal fragments)\nred = y ions (C-terminal fragments)\ngrey = noise",
        transform=ax.transAxes,ha="right",va="top",fontsize=9.4,color="#333",linespacing=1.5,
        bbox=dict(boxstyle="round,pad=0.42",fc="#FAFAFA",ec="#DDD"))
for i in range(len(bion)-1):
    ax.annotate("",xy=(bion[i+1],13),xytext=(bion[i],13),
                arrowprops=dict(arrowstyle="<|-|>",color="#7A8894",lw=0.9,shrinkA=1,shrinkB=1))
    ax.text((bion[i]+bion[i+1])/2,15.5,pep[i+1],ha="center",fontsize=8.6,color="#5A6570",fontweight="bold")
fig.text(0.5,-0.01,"The mass difference between consecutive fragments is one amino acid. That is how the sequence is read.",
        ha="center",fontsize=9.6,color="#555",style="italic")
fig.savefig(OUT+"bottomup_msms.png"); print("ok m2")
