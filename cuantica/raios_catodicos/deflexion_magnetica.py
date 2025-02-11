import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

U_as = 1e3 * np.vstack((1.1,1.5,2.0,2.5))
u_U_a = 1e3 * 0.1 / np.sqrt(12)

# 1.08e-2 ?
sq = 1e-2

sqx1 = np.array((9,8,6,5,4,3))
sqy1 = np.array((1,2,2,2,2,2))

sqx2 = np.array((9,8,7,7,6,5))
sqy2 = np.array((1,2,2,2,2,2))

sqx3 = np.array((8,7,6,5,4,3))
sqy3 = np.array((2,2,2,2,2,2))

sqx4 = np.array((8,7,6,5,4,3))
sqy4 = np.array((2,2,2,2,2,2))

sqxs = np.vstack((sqx1, sqx2, sqx3, sqx4))
sqys = np.vstack((sqy1, sqy2, sqy3, sqy4))

xs, ys = sqxs * sq, sqys * sq

rs = (xs**2 + ys**2) / (2*ys)

I1 = np.array((.07,.17,.27,.34,.45,.61))
I2 = np.array((.08,.19,.24,.31,.39,.53))
I3 = np.array((.23,.28,.36,.45,.60,.82))
I4 = np.array((.26,.32,.39,.52,.68,.93))

Is = np.vstack((I1,I2,I3,I4))
u_I = .01 / np.sqrt(12)

k = 4.2e-3

Bs = k * Is
u_B = k * u_I

#recta = lambda x, a, b : a + b*x
recta_origen = lambda x, b : b*x
cm = lambda b, U_a : 2*U_a/b
u_cm = lambda cm, b, U_a, u_b, u_U_a : cm * np.sqrt((u_b/b)**2 + (u_U_a/U_a)**2)

for i, U_a in enumerate(U_as):

    #popt, pcov = curve_fit(recta, 1/rs[i], Bs[i], sigma = u_B)
    popt, pcov = curve_fit(recta_origen, 1/rs[i]**2, Bs[i]**2, sigma = 2 * Bs[i] * u_B)
    
    #a, b = popt
    #u_a, u_b = np.sqrt(np.diag(pcov))

    b = popt[0]
    u_b = np.sqrt(np.diag(pcov))[0]

    pcm = cm(b, U_a)
    u_pcm = u_cm(pcm, b, U_a, u_b, u_U_a)

    #print(f"U_a = {U_a}\na = {a}({u_a})\nb = {b}({u_b})\ncm = {pcm}({u_pcm})\n-----")
    print(f"U_a = {U_a}\nb = {b}({u_b})\ncm = {pcm}({u_pcm})\n-----")

    fig, ax = plt.subplots()

    xlin = np.linspace(0.95*np.min(1/rs[i]**2), 1.05*np.max(1/rs[i]**2), 1000)

    #ax.plot(xlin, recta(xlin, a, b), ls = "solid", color = "tab:orange")
    ax.plot(xlin, 1e6 * recta_origen(xlin, b), ls = "solid", color = "tab:orange")
    ax.plot(1/rs[i]**2, 1e6 * Bs[i]**2, "o", color = "tab:blue") # (mT)**2

    ax.set_xlim(left = np.min(xlin), right = np.max(xlin))
    ax.set_ylim(bottom = 0)

    ax.set_xlabel(r"$r^{-2}$ (m$^{-2}$)")
    ax.set_ylabel(r"$B^2$ (mT$^2$)")

    fig.savefig(f"cuantica/raios_catodicos/deflexion_magnetica_{i}.pdf", dpi = 300, bbox_inches = "tight")
