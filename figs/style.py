import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

TEAL="#5F9EA0"; PURPLE="#9B8EC4"; OLIVE="#B5AE8A"; NAVY="#2E4057"
RED="#C1544B"; GREY="#B8B8B8"; LGREY="#E8E8E8"; ORANGE="#E08A3C"; BLUE="#4A77A8"

plt.rcParams.update({
    "font.family":"DejaVu Sans", "font.size":11,
    "axes.spines.top":False, "axes.spines.right":False,
    "axes.edgecolor":"#4A4A4A", "axes.labelcolor":"#2A2A2A",
    "axes.titlesize":13, "axes.titleweight":"bold", "axes.titlepad":12,
    "xtick.color":"#4A4A4A","ytick.color":"#4A4A4A",
    "figure.facecolor":"white","savefig.facecolor":"white",
    "savefig.dpi":160, "savefig.bbox":"tight",
})
OUT="/sessions/determined-friendly-gauss/mnt/lab-wiki/sources/07_textbooks/course_materials_2/assets/img/w7/"
