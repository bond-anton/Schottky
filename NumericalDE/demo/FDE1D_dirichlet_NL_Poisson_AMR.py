#-*- coding: utf-8 -*-
"""
Created on 30 мая 2015 г.

@author: anton
"""
from __future__ import division
import numpy as np


from matplotlib import pyplot as plt

from NumericalDE.FiniteDifference1D import dirichlet_non_linear_poisson_solver_amr

colors = ['b', 'g', 'y', 'k', 'm', 'c', 'b', 'g', 'y', 'k', 'm', 'c', 'b', 'g', 'y',
          'k', 'm', 'c', 'b', 'g', 'y', 'k', 'm', 'c']

Nd = lambda x: np.ones_like(x)
kT = 1/40

def f(x, Psi):    
    return 2*(1 - (np.exp(-Psi(x)/kT)))

def dfdDPsi(x, Psi):
    return 2/kT * np.exp(-Psi(x)/kT)

Psi = lambda x: np.exp(-0.7*x)
start = 0.0
stop = 20
bc1 = 1
bc2 = 0

N0 = 50

iteration = 1
idxs = np.array([0])
nodes = np.linspace(start, stop, num=int(N0*iteration) + 1, endpoint=True)

Meshes = dirichlet_non_linear_poisson_solver_amr(nodes, Psi, f, dfdDPsi, bc1, bc2,
                                                 max_iterations=1000, residual_threshold=1.5e-3, int_residual_threshold=1.5e-3,
                                                 max_level=20, mesh_refinement_threshold=1e-6, debug=True)

fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True)

flat_grid, flat_sol, flat_res = Meshes.flatten()

for level in Meshes.levels:
    for mesh in Meshes.Tree[level]:
        ax1.plot(mesh.phys_nodes(), mesh.solution, colors[level] + '-')
        dPsi = np.gradient(mesh.solution, mesh.phys_nodes(), edge_order=2)
        ax4.plot(mesh.phys_nodes(), dPsi, colors[level] + '-')

ax2.plot(flat_grid, flat_res, 'b-')
Meshes.plot_tree(ax3)
plt.show()
