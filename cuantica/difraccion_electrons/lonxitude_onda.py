import numpy as np

D1s = 1e-3 * np.array((29.845,27.650,25.870,24.915,22.550))
D2s = 1e-3 * np.array((51.280,46.800,43.185,40.305,38.455))

u_D_ind = 1e-3 * 0.1
u_D = u_D_ind / np.sqrt(2)

Vs = 1e3 * np.array((3.0,3.5,4.0,4.5,5.0))

u_V = 1e3 * 0.1 / np.sqrt(12)

L = 13.5e-2

d_1_th, d_2_th = 2.13e-10, 1.23e-10

lonxitude = lambda D, d : D*d / (2*L)
u_lonxitude = lambda lonx, D, u_D : lonx * u_D / D

lonx_1s = lonxitude(D1s, d_1_th)
u_lonx_1s = u_lonxitude(lonx_1s, D1s, u_D)
lonx_2s = lonxitude(D2s, d_2_th)
u_lonx_2s = u_lonxitude(lonx_2s, D2s, u_D)

q = 1.6021e-19
m = 9.1091e-31
h = 6.6256e-34

lonx_broglie = lambda V_A : h / np.sqrt(2*m*q*V_A)
u_lonx_broglie = lambda lonx_brog, V_A, u_V_A : lonx_brog * u_V_A / V_A / 2

lonx_broglies = lonx_broglie(Vs)
u_lonx_broglies = u_lonx_broglie(lonx_broglies, Vs, u_V)

for i, (l1, u_l1, l2, u_l2, l_b, u_l_b) in enumerate(zip(lonx_1s, u_lonx_1s, lonx_2s, u_lonx_2s, lonx_broglies, u_lonx_broglies)):
    print(f"V_A = {Vs[i]}")
    print(f"lonx_1 = {l1}({u_l1})")
    print(f"lonx_2 = {l2}({u_l2})")
    print(f"lonx_broglie = {l_b}({u_l_b})")
    print("-----")
