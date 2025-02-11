import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

U_as = 1e3 * np.vstack((1.0,1.5,2.0,2.5))
u_U_a = 1e3 * 0.1 / np.sqrt(12)

I1 = np.array((.15,.20,.30,.40,.50))
I2 = np.array((.13,.21,.33,.40,.54))
I3 = np.array((.15,.20,.30,.41,.50))
I4 = np.array((.15,.21,.33,.43,.51))

Is = np.vstack((I1,I2,I3,I4))
u_I = .01 / np.sqrt(12)

k = 4.2e-3
Bs = k * Is
u_B = k * u_I

U_p1 = 1e3 * np.array((0.6,0.8,1.0,1.2,1.3))
U_p2 = 1e3 * np.array((0.7,1.0,1.4,1.6,1.9))
U_p3 = 1e3 * np.array((0.9,1.1,1.5,1.9,2.2))
U_p4 = 1e3 * np.array((1.0,1.4,1.9,2.3,2.6))

U_ps = np.vstack((U_p1,U_p2,U_p3,U_p4))
u_U_p = 1e3 * 0.1 / np.sqrt(12)

d = 54e-3
Es = U_ps/ d
u_E = u_U_p / d

recta = lambda x, a, b : a + b*x
cm = lambda b, U_a : b**2 / 2 / U_a
u_cm = lambda cm, b, U_a, u_b, u_U_a : cm * np.sqrt((2*u_b/b)**2 + (u_U_a/U_a)**2)

for i, U_a in enumerate(U_as):

    popt, pcov = curve_fit(recta, Bs[i], Es[i], sigma = u_E)
    
    a, b = popt
    u_a, u_b = np.sqrt(np.diag(pcov))

    pcm = cm(b, U_a)
    u_pcm = u_cm(pcm, b, U_a, u_b, u_U_a)

    print(f"U_a = {U_a}\nv = {b}({u_b})\ncm = {pcm}({u_pcm})\n-----")

    fig, ax = plt.subplots()

    xlin = np.linspace(0.95*np.min(Bs[i]), 1.05*np.max(Bs[i]), 1000)

    ax.plot(xlin, recta(xlin, a, b), ls = "solid", color = "tab:orange")
    ax.plot(Bs[i], Es[i], "o", color = "tab:blue")

    ax.set_xlim(left = np.min(xlin), right = np.max(xlin))
    ax.set_ylim(bottom = 0)

    fig.savefig(f"cuantica/raios_catodicos/compensacion_campos_{i}.pdf", dpi = 300, bbox_inches = "tight")
