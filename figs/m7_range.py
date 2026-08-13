from style import *
import matplotlib.patches as mp
fig,ax=plt.subplots(figsize=(11.0,5.6))
prot=[("Albumin",10.6,"#C1544B"),("Immunoglobulins",10.0,"#C1544B"),
      ("Transferrin",9.2,"#E08A3C"),("Fibrinogen",9.0,"#E08A3C"),
      ("CRP",7.2,"#C9A227"),("PSA",4.6,TEAL),("Troponin I",3.2,TEAL),
      ("IL-6",2.4,BLUE),("TNF-α",2.0,BLUE),("Tumour-derived\nleakage proteins",1.0,PURPLE)]
for nm,v,c in prot:
    ax.barh([nm],[v],color=c,height=0.62,edgecolor="white")

ax.invert_yaxis()
ax.set_xlim(0,12.6)
ax.set_xlabel("plasma concentration,  log₁₀(pg / mL)")
ax.set_title("The plasma dynamic range problem",loc="left")
ax.axvspan(6.6,12.6,color="#C1544B",alpha=0.07,lw=0)
ax.text(9.6,-0.85,"top 20 proteins ≈ 99% of plasma protein mass",ha="center",fontsize=9.5,color="#8A3A33")
ax.annotate("",xy=(12.4,9.9),xytext=(6.7,9.9),arrowprops=dict(arrowstyle="<|-|>",color="#7A8894",lw=1.4))
ax.text(9.5,9.55,"unfractionated MS sees roughly this window (~5–6 orders)",
        ha="center",fontsize=9.3,color="#5A6570")
ax.annotate("",xy=(6.6,10.6),xytext=(0.3,10.6),arrowprops=dict(arrowstyle="<|-|>",color=TEAL,lw=1.6))
ax.text(3.4,10.25,"affinity assays (Olink, SomaScan) reach here",
        ha="center",fontsize=9.3,color="#2E5E4E",fontweight="bold")
ax.set_ylim(11.4,-1.4)
fig.text(0.5,-0.02,"This is why depleting albumin and immunoglobulins, or fractionating, is standard for plasma MS — and why affinity platforms are used for low-abundance biomarkers.",
         ha="center",fontsize=9.4,color="#666")
fig.savefig(OUT+"dynamic_range.png"); print("ok m7")
