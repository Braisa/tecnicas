import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

U_as = 1e3 * np.vstack((1.0,1.5,2.0,2.5))
u_U_a = 1e3 * 0.1 / np.sqrt(12)

U_p1 = np.array((0.5,0.7,0.8,1.0,1.2,1.6))
U_p2 = np.array((0.7,0.9,1.2,1.4,1.7,2.2))
U_p3 = np.array((0.9,1.2,1.5,1.8,2.2,2.8))
U_p4 = np.array((1.2,1.5,2.0,2.2,2.8,3.6))

U_ps = np.vstack((U_p1,U_p2,U_p3,U_p4))

# 1.08e-2 ?
sq = 1e-2

sqx = np.array((7,6,5,7,6,5))
sqy = np.array((1,1,1,2,2,2))

x = sq * sqx
y = sq * sqy

d = 54e-3

Es = U_ps / d

vs = np.sqrt(2*1.76e11*U_as)

# Non se fai por axustes ?
