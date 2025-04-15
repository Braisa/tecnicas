import numpy as np

D1s = 1e-3 * np.array((29.845,27.650,25.870,24.915,22.550))
D2s = 1e-3 * np.array((51.280,46.800,43.185,40.305,38.455))

u_D = 1e-3 * 0.35

Vs = 1e3 * np.array((3.0,3.5,4.0,4.5,5.0))

u_V = 1e3 * 0.1 / np.sqrt(12)

q = 1.6021e-19
m = 9.1091e-31
h = 6.6256e-34
L = 13.5e-2
u_L = 0.05e-2

d = lambda V, D : h / np.sqrt(8*m*q*V) / np.sin(.5*np.arctan(D/2/L))
u_d = lambda d, V, D, u_V, u_D : d*np.sqrt((u_V/2/V)**2 + (D**2 * u_D**2 + L**2 * u_L**2)/(D**2 + 4*L**2)**2 / np.tan(.5*np.arctan(D/2/L))**2)

d1s = d(Vs, D1s)
u_d1s = u_d(d1s, Vs, D1s, u_V, u_D)
d2s = d(Vs, D2s)
u_d2s = u_d(d2s, Vs, D2s, u_V, u_D)

d1 = np.mean(d1s)
u_d1 = np.sqrt(np.sum(u_d1s**2))/5
d2 = np.mean(d2s)
u_d2 = np.sqrt(np.sum(u_d2s**2))/5
