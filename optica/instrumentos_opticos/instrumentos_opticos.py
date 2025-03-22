import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import tabulate as tab

d_0 = 40
s_d_0 = 1

beta_exp = lambda y, Y : Y/y
s_beta_exp = lambda beta, y, Y, sy, sY : beta*np.sqrt((sY/Y)**2 + (sy/y)**2)

beta_lupa = lambda d, fp : 1 + (d_0 - d)/fp
s_beta_lupa = lambda d, fp, sd, sfp : 1/fp*np.sqrt(s_d_0**2 + sd**2 + (d-d_0)**2*sfp**2)

beta_micro = lambda d, fp, l, lp : lp/l * (1 - (d_0 - d)/fp)
s_beta_micro = lambda beta, d, fp, l, lp, sd, sfp, sl, slp : np.sqrt(s_beta_exp(beta, l, lp, sl, slp)**2 + s_beta_lupa(d,fp,sd,sfp)**2)

data_lentes = pd.read_excel("optica/instrumentos_opticos/lentes.xlsx")

recta = lambda x, a, b : a + b*x
focal = lambda a : 1/a
s_focal = lambda a, sa : sa/a**2

p1, sp1 = curve_fit(recta, 1/data_lentes["s"][0:3], 1/data_lentes["sp"][0:3])
p2, sp2 = curve_fit(recta, 1/data_lentes["s"][3:6], 1/data_lentes["sp"][3:6])
p3, sp3 = curve_fit(recta, 1/data_lentes["s"][6:8], 1/data_lentes["sp"][6:8])

a1, sa1 = p1[0], np.diag(np.sqrt(sp1))[0]
a2, sa2 = p2[0], np.diag(np.sqrt(sp2))[0]
a3, sa3 = p3[0], np.diag(np.sqrt(sp3))[0]

f1, sf1 = focal(a1), s_focal(a1, sa1)
f2, sf2 = focal(a2), s_focal(a2, sa2)
f3, sf3 = focal(a3), s_focal(a3, sa3)

data_lupa = pd.read_excel("optica/instrumentos_opticos/lupa.xlsx")

data_lupa["bexp"] = beta_exp(data_lupa["y"], data_lupa["Y"])
data_lupa["bth"] = beta_lupa(data_lupa["d"], f2)

"""print(tab.tabulate(data_lupa[["d","bth","bexp"]],
                   tablefmt = "latex_raw",
                   showindex = False
))
"""

data_lupa["s_bexp"] = s_beta_exp(data_lupa["bexp"], data_lupa["y"], data_lupa["Y"], 5, 5)
data_lupa["s_bth"] = s_beta_lupa(data_lupa["d"], f2, 1, sf2)

data_micro = pd.read_excel("optica/instrumentos_opticos/microscopio.xlsx")

s_lp = lambda lp, l, f, sl, sf : lp/np.abs(l+f) * np.sqrt((f*sl/l)**2 + (l*sf/f)**2)

"""print(tab.tabulate(data_micro[["d","L","Lp","y","Y"]],
                   tablefmt = "latex_raw",
                   showindex = False
))"""

f3, sf3 = 4.999, 0.058

slp = s_lp(data_micro["Lp"], data_micro["L"], f3, 0.1, sf3)

data_micro["bexp"] = beta_exp(data_micro["y"], data_micro["Y"])
data_micro["s_bexp"] = s_beta_exp(data_micro["bexp"], data_micro["y"], data_micro["Y"], 1, 1)

data_micro["bth"] = beta_micro(data_micro["d"], f2, data_micro["L"], data_micro["Lp"])
data_micro["s_bth"] = s_beta_micro(data_micro["bth"], data_micro["d"], f2, data_micro["L"], data_micro["Lp"], 0.1, sf2, 0.1, slp)

print(tab.tabulate(data_micro[["d", "bth", "bexp"]],
                   tablefmt = "latex_raw",
                   showindex = False
))

