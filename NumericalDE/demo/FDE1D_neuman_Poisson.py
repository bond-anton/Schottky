# -*- coding: utf-8 -*-

"""
Created on 30 мая 2015 г.

@author: anton
"""

import numpy as np
from matplotlib import pyplot as plt

from NumericalDE.FiniteDifference1D import neuman_poisson_solver


def f(x):
    return 2 * np.pi * np.sin(2 * np.pi * x) + 2


def Y(x):
    return -np.sin(2 * np.pi * x) / (2 * np.pi) + x ** 2 + x


nodes = np.linspace(0.0, 1.0, num=91, endpoint=True)
bc1 = 0
bc2 = 2
integral = np.trapz(f(nodes), nodes)
print integral, bc2 - bc1
print np.allclose(integral, bc2 - bc1)

Y_sol, R = neuman_poisson_solver(nodes, f, bc1, bc2, Psi0=0, debug=True)
dx = np.gradient(nodes)
dY_sol = np.gradient(Y_sol, dx, edge_order=2)
d2Y_sol = np.gradient(dY_sol, dx, edge_order=2)
fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
ax1.plot(nodes[2:-2], f(nodes[2:-2]), 'r-')
# ax1.plot(nodes, d2Y_sol, 'b-')
ax2.plot(nodes[2:-2], Y_sol[2:-2], 'b-')
ax2.plot(nodes[2:-2], Y(nodes[2:-2]), 'r-')
# ax1.plot(nodes, Y_sol-R, 'g-o')
# ax1.plot(nodes, Y_sol+R, 'g-o')
ax3.plot(nodes[2:-2], R[2:-2], 'g-o')
# ax2.plot(nodes, f(nodes), 'r-')
# ax2.plot(nodes, d2Y_sol, 'b-')
plt.show()
