import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import tabulate as tab
from sklearn.metrics import r2_score

angle = lambda deg, min : deg + min/60
s_angle = 1/60

rad = lambda ang : ang * np.pi/180
s_rad = s_angle * np.pi/180

b1, bm1 = 165, 26
b2, bm2 = 46, 28

beta1, beta2 = angle(b1, bm1), angle(b2, bm2)
alpha = np.abs(beta1 - beta2) / 2
s_alpha = s_angle / np.sqrt(2)

ref = angle(284,1)

dev = pd.read_excel("optica/dispersion_material/desviacion.xlsx")

dev["red"] = np.abs(angle(dev["r"], dev["rm"]) - ref)
dev["yellow"] = np.abs(angle(dev["y"], dev["ym"]) - ref)
dev["green"] = np.abs(angle(dev["g"], dev["gm"]) - ref)
dev["blue"] = np.abs(angle(dev["b"], dev["bm"]) - ref)
dev["purple"] = np.abs(angle(dev["p"], dev["pm"]) - ref)

print(tab.tabulate(dev[["red","yellow","green","blue","purple"]],
                   tablefmt = "latex_raw",
                   showindex = False)
)

s_dev = np.sqrt(2) * s_angle

s_mean = s_dev / np.sqrt(5)
dev["means"] = (np.mean(dev["red"]),np.mean(dev["yellow"]),np.mean(dev["green"]),np.mean(dev["blue"]),np.mean(dev["purple"]))

dev["deltas"] = rad(dev["means"])
s_delta = s_mean * np.pi/180

alpha = rad(alpha)
s_alpha *= np.pi/180

ind = lambda delta : np.sin((delta + alpha)/2)/np.sin(alpha/2)
s_ind = lambda delta, s_delta : 1/(2*np.sin(alpha/2))*np.sqrt((s_delta*np.cos((delta+alpha)/2))**2 + (s_alpha*np.sin(delta/2)/np.sin(alpha/2))**2)

dev["index"] = ind(dev["deltas"])
dev["s_index"] = s_ind(dev["deltas"], s_delta)

print(tab.tabulate(dev[["index", "s_index"]],
                   tablefmt = "latex_raw",
                   showindex = False)
)

recta = lambda x, a, b : a + b*x
cauchy = lambda x, a, b : a + b/(x**2)

dev["lam"] = (616e-9,589e-9,568e-9,515e-9,498e-9)
dev["x"] = 1/(dev["lam"]**2)

popt_recta, pcov_recta = curve_fit(recta, dev["x"], dev["index"], sigma = dev["s_index"])
popt_cauchy, pcov_cauchy = curve_fit(cauchy, dev["lam"], dev["index"], sigma = dev["s_index"])

fig, ax = plt.subplots(figsize=(4,4))

lamlin = np.linspace(485, 630, 1000)

ax.plot(lamlin, cauchy(lamlin * 1e-9, *popt_cauchy), ls = "dashed", color = "tab:blue", label = "Ajuste")
#ax.plot(dev["lam"] * 1e9, dev["index"], "o", color = "tab:orange", label = "Medidas")
ax.errorbar(dev["lam"] * 1e9, dev["index"], yerr = 2*dev["s_index"], fmt = ".", capsize = 3, color = "tab:orange", label = "Medidas")

ax.set_xlabel(r"$\lambda$ (nm)")
ax.set_ylabel(r"$n$")

ax.set_xlim(left = np.min(lamlin), right = np.max(lamlin))

ax.legend(loc = "best")
ax.xaxis.set_major_locator(plt.MultipleLocator(30))
ax.xaxis.set_minor_locator(plt.MultipleLocator(5))
ax.yaxis.set_major_locator(plt.MultipleLocator(0.01))
ax.yaxis.set_minor_locator(plt.MultipleLocator(0.0025))
fig.tight_layout()
fig.tight_layout()
fig.savefig("optica/dispersion_material/cauchy.pdf", dpi = 300, bbox_inches = "tight")

fig, ax = plt.subplots(figsize=(4,4))

xlin = 1/(np.linspace(485e-9, 630e-9, 1000))**2

ax.plot(xlin * 1e-12, recta(xlin, *popt_recta), ls = "dashed", color = "tab:blue", label = "Ajuste")
#ax.plot(dev["x"] * 1e-12, dev["index"], "o", color = "tab:orange", label = "Medidas")
ax.errorbar(dev["x"] * 1e-12, dev["index"], yerr = 2*dev["s_index"], fmt = ".", capsize = 3, color = "tab:orange", label = "Medidas")

ax.set_xlabel(r"$x$ ($\mu$m$^{-2}$)")
ax.set_ylabel(r"$n$")

ax.set_xlim(left = np.min(xlin*1e-12), right = np.max(xlin*1e-12))

ax.legend(loc = "best")
ax.xaxis.set_major_locator(plt.MultipleLocator(0.5))
ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))
ax.yaxis.set_major_locator(plt.MultipleLocator(0.01))
ax.yaxis.set_minor_locator(plt.MultipleLocator(0.0025))
fig.tight_layout()
fig.savefig("optica/dispersion_material/recta.pdf", dpi = 300, bbox_inches = "tight")

print(r2_score(dev["index"], recta(dev["x"], *popt_recta)))
