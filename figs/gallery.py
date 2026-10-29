# Look of the HaLab slide gallery (kkonoo.github.io/slide-gallery), course colour crimson:
# one dark emphasis per figure, tints of the same colour for steps, dark-grey text,
# thin grey lines, Pretendard with hierarchy by weight (Bold · SemiBold · Regular).
from style import *
import glob, os
from matplotlib import font_manager as fm
import matplotlib.patches as mp

for f in glob.glob(os.path.expanduser("~/.fonts/Pretendard-*.otf")):
    fm.fontManager.addfont(f)
plt.rcParams.update({"font.family":["Pretendard","DejaVu Sans"],
                     "text.color":"#262626","axes.labelcolor":"#262626"})

CRIMSON="#B0283C"; TEXT="#262626"; TEXT2="#666666"; LINE="#C9CED6"; FAINT="#F2F3F5"

def tint(a):
    """Crimson at opacity a (0–1) over white, as a solid colour."""
    c=np.array([0xB0,0x28,0x3C]); w=np.array([255,255,255])
    return "#%02X%02X%02X"%tuple((w+(c-w)*a).round().astype(int))

def card(ax,x,y,w,h,fc="white",ec=LINE,lw=0.9,z=2):
    ax.add_patch(mp.Rectangle((x-w/2,y-h/2),w,h,fc=fc,ec=ec,lw=lw,zorder=z))

def header(ax,eyebrow,title,x=0,dy=0,x1=1):
    """Gallery slide header in axes coords: crimson eyebrow, bold message title, thin rule.
    dy (points) lifts it, e.g. above panel titles; x1 stretches the rule across later panels."""
    from matplotlib.transforms import offset_copy
    kw=dict(xycoords="axes fraction",textcoords="offset points",ha="left",va="bottom")
    ax.annotate(eyebrow,(x,1),xytext=(0,34+dy),fontsize=10.5,color=CRIMSON,**kw)
    ax.annotate(title,(x,1),xytext=(0,10+dy),fontsize=16,fontweight="bold",color=TEXT,**kw)
    ax.plot([x,x1],[1,1],transform=offset_copy(ax.transAxes,ax.figure,y=dy,units="points"),
            color="#E3E6EB",lw=0.9,clip_on=False)
