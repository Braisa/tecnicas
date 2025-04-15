import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.cm import plasma
import cmasher as cmr
from scipy.optimize import curve_fit
import pandas as pd

mpl.rcParams["mathtext.fontset"] = "stix"
mpl.rcParams["font.family"] = "STIXGeneral"

# Imports

path = "cuantica/franck_hertz_hg/"
filenames = ("U2_0_4V.csv", "U2_0_7V.csv", "U2_1V.csv",
             "U2_1_5V.csv", "U2_2V.csv", "U2_2_5V.csv")

files = (pd.read_csv(path + filename, sep = "\s+", decimal = ",", names = ("U", "I", "i")) for filename in filenames)

sub_plasma = cmr.get_sub_cmap(plasma, .2, .9)
colors = sub_plasma(np.linspace(0, 1, np.size(filenames)))
labels = (r"$U_2 = 0.4$ V",r"$U_2 = 0.7$ V",r"$U_2 = 1.0$ V",
          r"$U_2 = 1.5$ V",r"$U_2 = 2.0$ V",r"$U_2 = 2.5$ V")

fig, ax = plt.subplots(figsize = (4,4))
U_max = 0

for _i, (file, c, l) in enumerate(zip(files, colors, labels)):
    ax.plot(file["U"], file["I"], marker = ",", color = c, label = l)
    U_max = np.max((U_max, np.max(file["U"])))

ax.set_xlabel(r"$U_1$ \left(V\right)")
ax.set_ylabel(r"$I$ \left(nA\right)")

ax.set_yscale("function", functions = (lambda x : x**(1/3), lambda x : x**3))

bounds = (.25,.55,.85,1.25,1.75,2.25,2.75)
norm = mpl.colors.BoundaryNorm(bounds, sub_plasma.N)

cbar = fig.colorbar(mpl.cm.ScalarMappable(norm = norm, cmap = sub_plasma),
                    ax = ax, orientation = "vertical", label = r"$U_2$ \left(V\right)")
cbar.ax.set_yticks((0.4,0.7,1.05,1.5,2.0,2.5))
cbar.ax.set_yticklabels((0.4,0.7,1.0,1.5,2.0,2.5))
cbar.ax.minorticks_off()

ax.set_xlim(left = 0, right = U_max)
ax.set_ylim(bottom = 0.007)

ax.xaxis.set_major_locator(plt.MultipleLocator(10))
ax.xaxis.set_minor_locator(plt.MultipleLocator(2))
ax.yaxis.set_major_locator(plt.AutoLocator())
ax.yaxis.set_minor_locator(plt.MultipleLocator(1))

fig.tight_layout()
fig.savefig(path + "curvas.pdf", dpi = 300, bbox_inches = "tight")
