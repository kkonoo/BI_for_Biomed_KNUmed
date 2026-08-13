from style import *
import matplotlib.patches as mp
fig,axes=plt.subplots(1,2,figsize=(12.4,4.8),gridspec_kw={"wspace":0.22,"width_ratios":[1.15,1]})

ax=axes[0]; ax.axis("off")
prots=[("Protein A",0.78,BLUE),("Protein B",0.50,TEAL),("Protein C",0.22,"#C9A227")]
peps=[("pep 1",0.06),("pep 2",0.22),("pep 3",0.38),("pep 4",0.54),("pep 5",0.70),("pep 6",0.86)]
for nm,y,c in prots:
    ax.add_patch(mp.FancyBboxPatch((0.03,y-0.048),0.20,0.096,
                 boxstyle="round,pad=0.006,rounding_size=0.02",fc=c,ec="none",alpha=0.85))
    ax.text(0.13,y,nm,ha="center",va="center",color="white",fontsize=10.5,fontweight="bold")
for nm,x in peps:
    ax.add_patch(mp.FancyBboxPatch((x-0.055,0.02),0.11,0.075,
                 boxstyle="round,pad=0.004,rounding_size=0.015",fc="#E4E8EC",ec="#9AA5B1",lw=1))
    ax.text(x,0.058,nm,ha="center",va="center",fontsize=9,color="#444")
links=[(0,0),(0,1),(0,2),(1,2),(1,3),(1,4),(2,4),(2,5)]
for pi,pj in links:
    y=prots[pi][1]; c=prots[pi][2]; x=peps[pj][1]
    shared=sum(1 for a,b in links if b==pj)>1
    ax.plot([0.23,x],[y,0.10],color=RED if shared else c,lw=2.0 if shared else 1.4,
            ls="--" if shared else "-",alpha=0.9 if shared else 0.55)
ax.text(0.62,0.955,"solid  = peptide unique to one protein",fontsize=9.4,color="#555",ha="center")
ax.text(0.62,0.90,"dashed red = peptide shared between proteins",fontsize=9.4,color=RED,ha="center")
ax.set_xlim(0,1); ax.set_ylim(0,1)
ax.set_title("Peptides map to proteins many-to-many",loc="left",fontsize=12.5)

ax=axes[1]; ax.axis("off")
ax.set_title("⚠️ The protein inference problem",loc="left",fontsize=12.5,color=RED)
blocks=[("If you only detect pep 3","Protein A and Protein B are both\nconsistent with the data. Neither is\nproven present.",RED),
        ("Protein group","Search engines report indistinguishable\nproteins as one 'protein group', with a\nrepresentative accession.",NAVY),
        ("Why it matters","Your 'differentially abundant protein'\nmay be a group of paralogues, and the\nreal one may not be the one named.",NAVY),
        ("What to check","How many unique peptides support it?\nOne shared peptide is not a protein\nidentification.",TEAL)]
y=0.90
for t,d,c in blocks:
    ax.text(0.0,y,t,fontsize=10.8,fontweight="bold",color=c,va="top")
    ax.text(0.0,y-0.075,d,fontsize=9.5,color="#555",va="top",linespacing=1.5)
    y-=0.245
ax.set_xlim(0,1); ax.set_ylim(-0.02,1)
fig.savefig(OUT+"protein_inference.png"); print("ok m4")
