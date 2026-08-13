from style import *
rng=np.random.default_rng(6)
fig,axes=plt.subplots(1,3,figsize=(13.5,4.4),gridspec_kw={"wspace":0.30})

# --- panel 1: library size ---
ax=axes[0]
ns=6; base=rng.lognormal(4.0,1.5,4000)
libs=np.array([1.0,0.95,1.05,2.9,3.1,2.85])
data=[np.log2(base*l+1) for l in libs]
bp=ax.boxplot(data,patch_artist=True,widths=0.6,showfliers=False)
for p,c in zip(bp["boxes"],[GREY]*3+[RED]*3):
    p.set_facecolor(c); p.set_alpha(0.6); p.set_edgecolor("#666")
for k in ["medians","whiskers","caps"]:
    for l in bp[k]: l.set_color("#555")
ax.set_xticklabels(["C1","C2","C3","T1","T2","T3"])
ax.set_ylabel("log₂(count + 1)")
ax.set_title("① Raw counts",loc="left",fontsize=12)
ax.text(0.5,0.05,"Treated samples were sequenced ~3× deeper.\nEvery gene looks 'up'.",transform=ax.transAxes,
        ha="center",fontsize=9.3,color="#8A3A33",linespacing=1.5,
        bbox=dict(boxstyle="round,pad=0.4",fc="#FDF3F2",ec="#E8C4C0"))

ax=axes[1]
data2=[np.log2(base*l/l+1) for l in libs]
data2=[d+rng.normal(0,0.04,len(d)) for d in data2]
bp=ax.boxplot(data2,patch_artist=True,widths=0.6,showfliers=False)
for p,c in zip(bp["boxes"],[GREY]*3+[TEAL]*3):
    p.set_facecolor(c); p.set_alpha(0.6); p.set_edgecolor("#666")
for k in ["medians","whiskers","caps"]:
    for l in bp[k]: l.set_color("#555")
ax.set_xticklabels(["C1","C2","C3","T1","T2","T3"])
ax.set_ylabel("log₂(normalised count + 1)")
ax.set_title("② After normalisation",loc="left",fontsize=12)
ax.text(0.5,0.05,"Distributions now comparable.\nThis is what CPM / TMM / RLE do.",transform=ax.transAxes,
        ha="center",fontsize=9.3,color="#2E5E4E",linespacing=1.5,
        bbox=dict(boxstyle="round,pad=0.4",fc="#F1F7F4",ec="#BFD8CC"))

# --- panel 3: composition effect ---
ax=axes[2]
genes=["gene A","gene B","gene C","gene D","GENE X\n(huge, treated only)"]
ctrl=np.array([100,100,100,100,10])
trt =np.array([100,100,100,100,900])
cs=ctrl/ctrl.sum()*1000; ts=trt/trt.sum()*1000
X=np.arange(5); w=0.36
ax.bar(X-w/2,cs,w,color=GREY,label="Control (as % of library)",edgecolor="white")
ax.bar(X+w/2,ts,w,color=RED,alpha=0.75,label="Treated (as % of library)",edgecolor="white")
ax.set_xticks(X); ax.set_xticklabels(genes,fontsize=8.6)
ax.set_ylabel("CPM")
ax.set_title("③ ⚠️ The composition effect",loc="left",fontsize=12,color=RED)
ax.legend(frameon=False,fontsize=8.8,loc="upper left")
ax.annotate("one gene explodes …",xy=(4,ts[4]),xytext=(3.0,ts[4]*0.86),fontsize=9,color="#8A3A33",
            ha="right",arrowprops=dict(arrowstyle="-|>",color="#C1544B",lw=1.1))
ax.annotate("… so every other gene\nfalls in relative terms,\nwith no change at all",
            xy=(1.18,ts[1]),xytext=(1.2,ts[4]*0.52),fontsize=9,color="#8A3A33",ha="center",
            linespacing=1.4,arrowprops=dict(arrowstyle="-|>",color="#C1544B",lw=1.1))
fig.text(0.5,-0.04,"CPM corrects depth only. TMM (edgeR) and RLE (DESeq2) also correct composition — which is why they are the standard for differential expression.",
         ha="center",fontsize=9.3,color="#666")
fig.savefig(OUT+"normalization.png"); print("ok h2")
