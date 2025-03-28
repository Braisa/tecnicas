import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

D1s = 1e-3 * np.array((29.845,27.650,25.870,24.915,22.550))
D2s = 1e-3 * np.array((51.280,46.800,43.185,40.305,38.455))

u_D_ind = 1e-3 * 0.1
u_D = u_D_ind / np.sqrt(2)

Vs = 1e3 * np.array((3.0,3.5,4.0,4.5,5.0))

u_V = 1e3 * 0.1 / np.sqrt(12)

recta = lambda x, k : k*x

k1, u_k1 = curve_fit(recta, 1/np.sqrt(Vs), D1s, sigma = u_D)
k2, u_k2 = curve_fit(recta, 1/np.sqrt(Vs), D2s, sigma = u_D)

fig, ax = plt.subplots()

xlin = np.linspace(.95*np.min(1/np.sqrt(Vs)), 1.05*np.max(1/np.sqrt(Vs)), 1000)

ax.plot(xlin, 1e3 * recta(xlin, k1), ls = "dashed", color = "tab:orange", label = " ")
ax.plot(xlin, 1e3 * recta(xlin, k2), ls = "dashed", color = "tab:green", label = " ")
ax.plot(1/np.sqrt(Vs), 1e3 * D1s, "o", color = "tab:orange", label = r"$k_1$")
ax.plot(1/np.sqrt(Vs), 1e3 * D2s, "^", color = "tab:green", label = r"$k_2$")

ax.set_xlim(left = np.min(xlin), right = np.max(xlin))
ax.legend(loc = "best", ncol = 2, handletextpad = 0.3, columnspacing = -0.6)

ax.set_xlabel(r"$V_A^{-\dfrac{1}{2}}$ (V$^{-\dfrac{1}{2}}$)")
ax.set_ylabel(r"$D$ (mm)")

fig.savefig("cuantica/difraccion_electrons/determinacion_lattice.pdf", dpi = 300, bbox_inches = "tight")

print(f"k1 = {k1}({u_k1})")
print(f"k2 = {k2}({u_k2})")

q = 1.6021e-19
m = 9.1091e-31
h = 6.6256e-34
L = 13.5e-2

d = lambda k : 2*L*h/(k*np.sqrt(2*m*q))
u_d = lambda d, k , u_k : d * u_k / k

d1, d2 = d(k1), d(k2)
u_d1, u_d2 = u_d(d1, k1, u_k1), u_d(d2, k2, u_k2)

print(f"d1 = {d1}({u_d1})")
print(f"d2 = {d2}({u_d2})")
