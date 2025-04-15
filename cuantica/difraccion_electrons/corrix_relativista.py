import numpy as np

Vs = 1e3 * np.array((3.0,3.5,4.0,4.5,5.0))
u_V = 1e3 * 0.1 / np.sqrt(12)

q = 1.6021e-19
m = 9.1091e-31
c = 299_792_458

corrix = lambda V : 1 - (q*V/m/c**2)
u_corrix = lambda corrix, V : corrix*u_V/V

alphas = corrix(Vs)
u_alphas = u_corrix(alphas, Vs)

l1s = 1e-12 * np.array((23.54,21.81,20.41,19.66,17.79))
u_l1s = 1e-14 * np.array((29,29,29,29,29))

l2s = 1e-12 * np.array((23.36,21.32,19.67,18.36,17.52))
u_l2s = 1e-14 * np.array((18,18,18,17,17))

l_corrix = lambda l, alpha : l*alpha
u_l_corrix = lambda l_corrix, l, alpha, u_l, u_alpha : l_corrix*np.sqrt((u_l/l)**2 + (u_alpha/alpha)**2)

l1s_corrix = l_corrix(l1s, alphas)
u_l1s_corrix = u_l_corrix(l1s_corrix, l1s, alphas, u_l1s, u_alphas)

l2s_corrix = l_corrix(l2s, alphas)
u_l2s_corrix = u_l_corrix(l2s_corrix, l2s, alphas, u_l2s, u_alphas)
