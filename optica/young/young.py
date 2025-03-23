import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import tabulate as tab
from sklearn.metrics import r2_score

yg = pd.read_excel("optica/young/young.xlsx")

d_sensor = .32
s_D = (np.sqrt(2) * 0.1) * 1e-2

yg["Dr1"] = (yg["D1"] - d_sensor) * 1e-2
yg["Dr2"] = (yg["D2"] - d_sensor) * 1e-2

yg["Deltar1"] = yg["Delta1"] * 5.2e-6
yg["Deltar2"] = yg["Delta2"] * 5.2e-6

s_Deltar = 10 * 5.2e-6

yg["i1"] = yg["Deltar1"]/(yg["N1"]-1)
yg["i2"] = yg["Deltar2"]/(yg["N2"]-1)

yg["s_i1"] = s_Deltar/(yg["N1"]-1)
yg["s_i2"] = s_Deltar/(yg["N2"]-1)

d1s = (72.01, 72.00, 73.00)
d2s = (36.01, 36.00, 37.00)
s_ds = .01

d1, d2 = np.mean(d1s), np.mean(d2s)
s_d = s_ds/np.sqrt(3)

beta, s_beta = 0.729, 0#.087
d1r = d1 * 5.2e-6 / beta
d2r = d2 *5.2e-6 / beta
s_d1r = d1r * np.sqrt((s_beta/beta)**2 + (s_d/d1)**2)
s_d2r = d2r * np.sqrt((s_beta/beta)**2 + (s_d/d2)**2)

yg["lam1"] = yg["i1"] * d1r / yg["Dr1"]
yg["lam2"] = yg["i2"] * d2r / yg["Dr2"]

yg["s_lam1"] = yg["lam1"] * np.sqrt((yg["s_i1"]/yg["i1"])**2 + (s_d1r/d1r)**2 + (s_D/yg["Dr1"])**2)
yg["s_lam2"] = yg["lam2"] * np.sqrt((yg["s_i2"]/yg["i2"])**2 + (s_d2r/d2r)**2 + (s_D/yg["Dr2"])**2)

recta = lambda x, b : b*x
lam = lambda b, d : b*d
s_lam = lambda lam, b, d, sb, sd : lam*np.sqrt((sb/b)**2 + (sd/d)**2)

b1, c1 = curve_fit(recta, yg["Dr1"].drop(3), yg["i1"].drop(3), sigma = yg["s_i1"].drop(3))
b2, c2 = curve_fit(recta, yg["Dr2"].drop(1), yg["i2"].drop(1), sigma = yg["s_i2"].drop(1))

sb1, sb2 = np.diag(np.sqrt(c1)), np.diag(np.sqrt(c2))

lam1, lam2 = lam(b1,d1r), lam(b2,d2r)
s_lam1, s_lam2 = s_lam(lam1,b1,d1r,sb1,s_d1r), s_lam(lam2,b2,d2r,sb2,s_d2r)

r1 = r2_score(yg["i1"].drop(3), recta(yg["Dr1"].drop(3), b1))
r2 = r2_score(yg["i2"].drop(1), recta(yg["Dr2"].drop(1), b2))

fig, ax = plt.subplots(figsize=(6,12))

Dlin = np.linspace(30.5, 48.5, 1000)

ax.plot(Dlin, recta(Dlin, b1*10), ls = "dashed", color = "tab:blue", label = "Ajuste Pareja I")
ax.plot(Dlin, recta(Dlin, b2*10), ls = "dashed", color = "tab:purple", label = "Ajuste Pareja II")

ax.errorbar(yg["Dr1"].drop(3) * 1e2, yg["i1"].drop(3) * 1e3, yerr = yg["s_i1"].drop(3) * 1e3, fmt = ".", capsize = 3, color = "tab:orange", label = "Medidas Pareja I")
ax.errorbar(yg["Dr2"].drop(1) * 1e2, yg["i2"].drop(1) * 1e3, yerr = yg["s_i2"].drop(1) * 1e3, fmt = ".", capsize = 3, color = "tab:pink", label = "Medidas Pareja II")

ax.set_xlabel(r"$D$ (cm)")
ax.set_ylabel(r"$i$ (mm)")

ax.set_xlim(left = np.min(Dlin), right = np.max(Dlin))

ax.legend(loc = "best")

plt.tight_layout()
fig.savefig("optica/young/young.pdf", dpi = 300, bbox_inches = "tight")
