import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.cm import plasma
import cmasher as cmr
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
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
labels_readable = ("0_4","0_7","1_0","1_5","2_0","2_5")

fig, ax = plt.subplots(figsize = (4,4))
fig_v, ax_v = plt.subplots(figsize = (4,4))
U_max = 0

# For peaks

parab = lambda x, a, b, c : a + b*x + c*x**2
vertex = lambda b, c : -b/2/c
u_vertex = lambda vertix, b, c, u_b, u_c : vertix * np.sqrt((u_b/b)**2 + (u_c/c)**2)

peak_min_heights = (2.2,2,1.2,0.4,0.2, 0.04)
distances = (50,50,75,70,75,75)
peak_half_width = 25

# For peak distances
u_distance = lambda u1, u2 : np.sqrt(u1**2 + u2**2)

for i, (file, cl, l, lr) in enumerate(zip(files, colors, labels, labels_readable)):
    ax.plot(file["U"], file["I"], marker = ",", color = cl, label = l)
    ax_v.plot(file["U"], file["I"], marker = ",", color = cl, label = l)
    U_max = np.max((U_max, np.max(file["U"])))
    
    peaks, _ = find_peaks(file["I"], height = peak_min_heights[i], distance = distances[i])
    #ax.plot(file["U"][peaks], file["I"][peaks], "o", color = cl)

    deltas_maxs, u_deltas_maxs = np.array(()), np.array(())
    deltas_mins, u_deltas_mins = np.array(()), np.array(())

    print(f"Curva {l}\n--------")
    for j, peak in enumerate(peaks):
        neighboring_pos = peak + np.arange(-peak_half_width, peak_half_width)
        popt, pcov = curve_fit(parab, file["U"][neighboring_pos], file["I"][neighboring_pos])
        
        a, b, c = popt
        u_a, u_b, u_c = np.diag(np.sqrt(pcov))

        vert = vertex(b, c)
        u_vert = u_vertex(vert, b, c, u_b, u_c)

        match j % 2:
            case 0:
                name = "max" + str(int((j+2)/2))
                if j > 0:
                    deltas_maxs = np.append(deltas_maxs, vert-prev_max)
                    u_deltas_maxs = np.append(u_deltas_maxs, u_distance(u_prev_max, u_vert))
                prev_max, u_prev_max = vert, u_vert
            case 1:
                name = "min" + str(int(1+(j-1)/2))
                if j > 1:
                    deltas_mins = np.append(deltas_mins, vert-prev_min)
                    u_deltas_mins = np.append(u_deltas_mins, u_distance(u_prev_min, u_vert))
                prev_min, u_prev_min = vert, u_vert
        """
        print(f"Peak {j}")
        print(f"a = {a}({u_a})")
        print(f"b = {b}({u_b})")
        print(f"c = {c}({u_c})")
        print(f"v = {vert}({u_vert})")
        print("\n")
        """
        #ax_v.plot(file["U"][neighboring_pos], parab(file["U"][neighboring_pos], *popt), color = "tab:gray")
        ax_v.plot(vert, parab(vert, *popt), "o", color = "tab:gray")

        fig_peak, ax_peak = plt.subplots(figsize = (4,4))

        ax_peak.plot(file["U"][neighboring_pos], file["I"][neighboring_pos], "x", color = cl, label = "Medidas")
        ax_peak.plot(file["U"][neighboring_pos], parab(file["U"][neighboring_pos], *popt), ls = "solid", color = cl, label = "Axuste")
        ax_peak.plot(vert, parab(vert, *popt), "o", color = cl, label = "Vértice")

        ax_peak.set_xlim(left = .99*np.min(file["U"][neighboring_pos]), right = 1.01*np.max(file["U"][neighboring_pos]))
        if lr == "2_5" and name == "min3":
            ax_peak.set_xlim(right = 1.01*vert)

        ax_peak.set_xlabel(r"$U_1$ (V)")
        ax_peak.set_ylabel(r"$I$ (nA)")

        ax_peak.legend(loc = "best")

        fig_peak.tight_layout()
        fig_peak.savefig(path + f"curva_{lr}/{lr}_" + name + ".pdf", dpi = 300, bbox_inches = "tight")

    """
    max_md = np.average(deltas_maxs, weights = u_deltas_maxs**-2)
    u_max_md = 1/np.sqrt(np.sum(u_deltas_maxs**-2))
    min_md = np.average(deltas_mins, weights = u_deltas_mins**-2)
    u_min_md = 1/np.sqrt(np.sum(u_deltas_mins**-2))
    print(f"max_md = {max_md}({u_max_md})")
    print(f"min_md = {min_md}({u_min_md})")
    """

ax.set_xlabel(r"$U_1$ (V)")
ax.set_ylabel(r"$I$ (nA)")

ax.set_yscale("function", functions = (lambda x : x**(1/3), lambda x : x**3))

bounds = (.25,.55,.85,1.25,1.75,2.25,2.75)
norm = mpl.colors.BoundaryNorm(bounds, sub_plasma.N)

cbar = fig.colorbar(mpl.cm.ScalarMappable(norm = norm, cmap = sub_plasma),
                    ax = ax, orientation = "vertical", label = r"$U_2$ (V)")
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

ax_v.set_xlabel(r"$U_1$ (V)")
ax_v.set_ylabel(r"$I$ (nA)")

ax_v.set_yscale("function", functions = (lambda x : x**(1/3), lambda x : x**3))

cbar = fig.colorbar(mpl.cm.ScalarMappable(norm = norm, cmap = sub_plasma),
                    ax = ax_v, orientation = "vertical", label = r"$U_2$ (V)")
cbar.ax.set_yticks((0.4,0.7,1.05,1.5,2.0,2.5))
cbar.ax.set_yticklabels((0.4,0.7,1.0,1.5,2.0,2.5))
cbar.ax.minorticks_off()

ax_v.set_xlim(left = 0, right = U_max)
ax_v.set_ylim(bottom = 0.007)

ax_v.xaxis.set_major_locator(plt.MultipleLocator(10))
ax_v.xaxis.set_minor_locator(plt.MultipleLocator(2))
ax_v.yaxis.set_major_locator(plt.AutoLocator())
ax_v.yaxis.set_minor_locator(plt.MultipleLocator(1))

fig_v.tight_layout()
fig_v.savefig(path + "curvas_vertex.pdf", dpi = 300, bbox_inches = "tight")
