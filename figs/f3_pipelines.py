from gallery import *

fig, ax = plt.subplots(figsize=(11.5,7.4))
ax.set_xlim(0,1); ax.set_ylim(0.04,0.95); ax.axis("off")
header(ax,"From raw data to variants","Two pipelines, one destination: a VCF")

lanes=[(0.0,"Sequencing: WGS / WES","GATK Best Practices — starts from reads",
        [("Sequencer output","reads + per-base quality","FASTQ"),
         ("Align to the reference","BWA-MEM, then sort","BAM"),
         ("Mark duplicates","flag PCR copies of one fragment",None),
         ("Recalibrate base qualities","BQSR",None),
         ("Call variants per sample","HaplotypeCaller, GVCF mode","gVCF"),
         ("Joint genotyping","GenomicsDBImport → GenotypeGVCFs",None),
         ("Filter","VQSR or hard filters","VCF")]),
       (0.52,"Genotyping array","no reads — starts from fluorescence",
        [("Scanner output","two-colour intensity per probe","IDAT"),
         ("Genotype calling","cluster intensities → AA / AB / BB","PLINK"),
         ("Sample QC","call rate · sex · heterozygosity · relatedness",None),
         ("Variant QC","call rate · HWE · MAF",None),
         ("Match the reference panel","strand · genome build · allele order",None),
         ("Phasing","SHAPEIT / Eagle",None),
         ("Imputation","Minimac / Beagle  (§2.4)","VCF + dosage")])]

LW,top,dy=0.48,0.80,0.085
for x0,title,sub,steps in lanes:
    ax.text(x0,0.905,title,ha="left",va="center",fontsize=13,fontweight="bold",color=TEXT)
    ax.text(x0,0.872,sub,ha="left",va="center",fontsize=9.5,color=TEXT2)
    cx=x0+0.02; ylast=top-(len(steps)-1)*dy
    ax.plot([cx,cx],[top,ylast],color=LINE,lw=1.1,ls=(0,(2,2.5)),zorder=1)
    for i,(name,tool,out) in enumerate(steps):
        y=top-i*dy
        ax.scatter([cx],[y],s=430,color=tint(0.16),zorder=2)
        ax.text(cx,y,str(i+1),ha="center",va="center",fontsize=10,fontweight="bold",color=CRIMSON,zorder=3)
        ax.text(cx+0.03,y+0.014,name,ha="left",va="center",fontsize=11,fontweight="semibold",color=TEXT)
        ax.text(cx+0.03,y-0.017,tool,ha="left",va="center",fontsize=9.2,color=TEXT2)
        if out:
            ax.text(x0+LW-0.01,y,out,ha="right",va="center",fontsize=9,color=TEXT,family="DejaVu Sans Mono",
                    bbox=dict(boxstyle="round,pad=0.35,rounding_size=0.8",fc="white",ec=LINE,lw=0.9))
    # into the shared output
    ax.plot([cx,cx,0.5],[ylast-0.035,0.175,0.175],color="#9AA1AB",lw=1,zorder=1)
ax.plot([0.5,0.5],[0.175,0.15],color="#9AA1AB",lw=1,zorder=1)
ax.plot([0.5,0.5],[0.85,0.27],color="#E3E6EB",lw=0.9)        # lane divider

card(ax,0.5,0.105,0.52,0.085,fc=CRIMSON,ec=CRIMSON)
ax.text(0.5,0.118,"VCF: a genotype for every sample at every site",ha="center",va="center",
        color="white",fontsize=12,fontweight="bold",zorder=3)
ax.text(0.5,0.086,"→ association (§3)  ·  clinical interpretation (§6)",ha="center",va="center",
        color=tint(0.18),fontsize=9.5,zorder=3)
fig.savefig(OUT+"pipelines.png"); print("ok f3 pipelines")
