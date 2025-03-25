import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score

ns = np.array((-3.,-2.,-1.,1.,2.,3.), dtype = float)
lam = 632.8e-9

s_x = 0.1

x11 = 1e-1 * np.array((-85,-50,-26,24,48,73))
x12 = 1e-1 * np.array((-84,-54,-27,28,56,83))
x13 = 1e-1 * np.array((-75,-49,-25,24,49,73))
x14 = 1e-1 * np.array((-62,-40,-20,21,41,60))

x21 = 1e-1 * np.array((-49,-22,-17,17,33,49))
x22 = 1e-1 * np.array((-46,-31,-17,14,29,43))
x23 = 1e-1 * np.array((-42,-23,-14,14,24,42))
x24 = 1e-1 * np.array((-37,-25,-13,14,26,37))

s_z = 0.5
zs = (304.5,342.0,358.0,390.0,241.0,290.5,324.0,356.0)

recta = lambda n, x0, b : x0 + b*n
ancho = lambda b, z : lam*z/b
s_ancho = lambda a, b, z, sb, sz : a*np.sqrt((sz/z)**2 + (sb/b)**2)

xs = (x11,x12,x13,x14,x21,x22,x23,x24)

labels = ("z11","z12","z13","z14","z21","z22","z23","z24")

for i, (z, x, l) in enumerate(zip(zs,xs,labels)):

    popt, pcov = curve_fit(recta, ns, x, sigma = s_x)

    x0, b = popt
    s_x0, s_b = np.diag(np.sqrt(pcov))

    r2 = r2_score(x, recta(ns, x0, b))

    a = ancho(b, z)
    s_a = s_ancho(a, b, z, s_b, s_z)

    print(f"{l}")
    print(f"a = {a}({s_a})")
    print(f"x0 = {x0}({s_x0})")
    print(f"b = {b}({s_b})")
    print(f"r2 = {r2}")
    print("\n")

    fig, ax = plt.subplots()

    nlin = np.linspace(-3.1, 3.1, 100)

    ax.plot(nlin, recta(nlin, *popt), ls = "dashed", color = "tab:blue", label = "Ajuste")
    ax.plot(ns, x, "o", color = "tab:orange", label = "Medidas")
    
    ax.set_xlim(left = -3.1, right = 3.1)
    
    ax.set_xlabel(r"$n$")
    ax.set_ylabel(r"$x_\text{mín}$ (cm)")

    ax.legend(loc = "best")
    fig.tight_layout()
    fig.savefig("optica/difraccion/" + l + ".pdf", dpi = 300, bbox_inches = "tight")

