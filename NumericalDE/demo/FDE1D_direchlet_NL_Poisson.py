# -*- coding: utf-8 -*-
"""
Created on 06 июня 2015 г.

@author: anton
"""
from __future__ import division
import numpy as np
from matplotlib import pyplot as plt

from NumericalDE.FiniteDifference1D import dirichlet_non_linear_poisson_solver

Nd = lambda x: np.ones_like(x)
kT = 1 / 20


def f(x, Psi):
    return Nd(x) * (1 - (np.exp(-Psi(x) / kT)))


def dfdDPsi(x, Psi):
    return Nd(x) / kT * np.exp(-Psi(x) / kT)


Psi = lambda x: np.exp(-x * 3)

nodes = np.linspace(0., 4., num=51, endpoint=True)
bc1 = 1
bc2 = 0

DPsi = np.zeros_like(nodes)
E = np.zeros_like(nodes)
plt.ion()
_, (ax1, ax2, ax3, ax4) = plt.subplots(4)
ax1.set_autoscaley_on(True)
ax2.set_autoscaley_on(True)
ax3.set_autoscaley_on(True)
ax4.set_autoscaley_on(True)
Psi_line, = ax1.plot(nodes, Psi(nodes))
DPsi_line, = ax2.plot(nodes, DPsi)
dPsi = np.gradient(Psi(nodes), nodes, edge_order=2)
d2Psi = np.gradient(dPsi, nodes, edge_order=2)
f_line, = ax3.plot(nodes, f(nodes, Psi))
d2Psi_line, = ax3.plot(nodes, d2Psi)
E_line, = ax4.plot(nodes, E)

plt.draw()
# plt.show()

for i in range(1000):
    print i + 1
    # time.sleep(0.1)
    Psi, DPsi, R = dirichlet_non_linear_poisson_solver(nodes, Psi, f, dfdDPsi, bc1=1, bc2=0, J=1, debug=False)
    dPsi = np.gradient(Psi(nodes), nodes, edge_order=2)
    d2Psi = np.gradient(dPsi, nodes, edge_order=2)
    Psi_line.set_ydata(Psi(nodes))
    DPsi_line.set_ydata(DPsi)
    f_line.set_ydata(f(nodes, Psi))
    d2Psi_line.set_ydata(d2Psi)
    E_line.set_ydata(R)
    ax1.relim()
    ax1.autoscale_view()
    ax2.relim()
    ax2.autoscale_view()
    ax3.relim()
    ax3.autoscale_view()
    ax4.relim()
    ax4.autoscale_view()
    plt.draw()
