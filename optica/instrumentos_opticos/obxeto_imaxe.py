import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

s = -1e-2 * np.array((15,15,15,15,
                      20,20,20,20,
                      18,18,18,18,
                      13,13,13,13))

u_s = 1e-2 * 0.1

sp = 1e-2 * np.array((32.9,33.3,33.1,33.5,
                      21.3,21.3,21.4,21.1,
                      24.0,24.1,24.3,24.0,
                      49.5,50.0,49.5,50.3))

u_sp = 1e-2 * 0.1

recta = lambda x, a, b : a + b*x

x = 1/s
u_x = u_s / s**2

popt, pcov = curve_fit(recta, x, 1/sp, sigma = u_x)

fp, slope = 1/popt[0], popt[1]
u_fp, u_slope = fp**2 * np.sqrt(np.diag(pcov))[0], np.sqrt(np.diag(pcov))[1]

print(f"fp = {fp}({u_fp})")
print(f"slope = {slope}({u_slope})")

fig, ax = plt.subplots()

x_lin = -1 * np.linspace(4.5, 8, 1000)

ax.plot(x_lin, recta(x_lin, *popt), ls = "solid", color = "tab:orange")
ax.plot(x, 1/sp, "o", color = "tab:blue")

ax.set_xlim(left = np.min(x_lin), right = np.max(x_lin))
ax.set_ylim(bottom = 0.95 * np.min(1/sp), top = 1.05 * np.max(1/sp))

ax.set_xlabel(r"$s^{-1}$ (m$^{-1}$)")
ax.set_ylabel(r"$s\prime^{-1}$ (m$^{-1}$)")

fig.savefig("optica/instrumentos_opticos/obxeto_imaxe.pdf", dpi = 300, bbox_inches = "tight")
