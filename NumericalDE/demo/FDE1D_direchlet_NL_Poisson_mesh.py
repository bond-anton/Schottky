# -*- coding: utf-8 -*-
"""
Created on 06 июня 2015 г.

@author: anton
"""
from __future__ import division
import numpy as np
from matplotlib import pyplot as plt

from NumericalDE.FiniteDifference1D import dirichlet_non_linear_poisson_solver_reccurent_mesh
from NumericalDE.Mesh import UniformMesh1D, Uniform1DMeshesTree

Nd = lambda x: np.ones_like(x)
kT = 1 / 20


def f(x, Psi):
    return Nd(x) * (1 - (np.exp(-Psi(x) / kT)))


def dfdDPsi(x, Psi):
    return Nd(x) / kT * np.exp(-Psi(x) / kT)


Psi = lambda x: np.exp(-x * 3)

nodes = np.linspace(0., 10., num=51, endpoint=True)
bc1 = 1
bc2 = 0

root_mesh = UniformMesh1D(nodes[0], nodes[-1], nodes[1] - nodes[0], bc1, bc2)

Meshes = Uniform1DMeshesTree(root_mesh, refinement_coeficient=2, aligned=True)

root_mesh, Psi = dirichlet_non_linear_poisson_solver_reccurent_mesh(root_mesh, Psi, f, dfdDPsi, max_iterations=1000,
                                                                    threshold=1e-6, debug=True)

mesh_refinement_threshold = 1e-7

idxs = np.where(abs(root_mesh.residual) > mesh_refinement_threshold)

dx = np.gradient(root_mesh.phys_nodes())
dPsi = np.gradient(root_mesh.solution, dx, edge_order=2)
d2Psi = np.gradient(dPsi, dx, edge_order=2)

_, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
ax1.plot(root_mesh.phys_nodes(), root_mesh.solution)
ax1.plot(root_mesh.phys_nodes()[idxs], root_mesh.solution[idxs], 'radius-o')
ax2.plot(root_mesh.phys_nodes(), root_mesh.residual)
ax2.plot(root_mesh.phys_nodes()[idxs], root_mesh.residual[idxs], 'radius-o')
ax3.plot(root_mesh.phys_nodes(), f(root_mesh.phys_nodes(), Psi))
ax3.plot(root_mesh.phys_nodes(), d2Psi)
plt.show()
