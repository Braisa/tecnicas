import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pch

theta = np.radians(10 * np.arange(36))
P = 1e-2 * np.array((81,86,90,89,81,82,73,64,52,43,35,31,31,34,40,49,58,68,
                     82,92,97,94,94,85,75,62,50,41,33,29,28,33,41,50,60,70))
a1, b1, o1 = .95, .6, 25

fig, ax = plt.subplots()

theta_lin = np.radians(np.linspace(0, 360, 100))

ax.add_patch(pch.Ellipse((0,0),2*a1,2*b1,angle = o1, fill = False, color = "tab:purple", label = r"$\theta_0 = 25^\circ$"))
ax.add_patch(pch.Ellipse((0,0),2*a1,2*b1,angle = 30, fill = False, color = "tab:blue", label = r"$\theta_0 = 30^\circ$"))

ax.plot(np.sqrt(P)*np.cos(theta), np.sqrt(P)*np.sin(theta), "o", color = "tab:orange", label = "Medidas")

ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$y$")

ax.legend(loc = "best")

fig.savefig("optica/polarizacion/elipse.pdf", dpi = 300, bbox_inches = "tight")
