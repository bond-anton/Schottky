# -*- coding: utf-8 -*-

'''
Created on 10 апр. 2015 г.

@author: anton
'''

from __future__ import division, print_function
from os.path import dirname, join

from pathos.pools import ProcessPool as Pool

import numpy as np

np.seterr(all='raise')
import mpmath as mp
# from scipy.interpolate import interp1d
# mp.mp.dps = 25

import matplotlib
# matplotlib.use("Qt4Agg")
from matplotlib import pyplot as plt

# plt.xkcd()

from ProjectManager import Project

from Schottky.Notation import q, Eg
from Schottky.Metal import Metal
from Schottky.Semiconductor_old import Semiconductor, Trap, Dopant, Dislocation, BondingInterface
from Schottky.Diode import SchottkyDiode, Poisson, Kinetics, Visual
from Schottky.Helpers import Psi_approx
from .diode_description import *

colors = ['b', 'g', 'y', 'k', 'm', 'c', 'b', 'g', 'y', 'k', 'm', 'c', 'b', 'g', 'y', 'k', 'm', 'c', 'b', 'g', 'y', 'k',
          'm', 'c']

'''
Electrode = Metal('Au', q * 5.1)

P = Dopant('Phosphorus', 1.0e21, [[+1, 0.045 * q, 1], [0, 0.045 * q, 1]])

Si = Semiconductor('Si', lookup=True)
Si.add_dopant(P)

a = 5.43e-10
b = a * np.sqrt(2) / 2

S_e1 = Trap('Shallow electron trap', charge_states=[[0, 0.1 * q, 1], [-1, 0.1 * q, 1]],
            energy_distribution_function='Single Level', energy_spread=0.01 * q)
D_e1 = Trap('Deep electron trap', charge_states=[[0, 0.3 * q, 1], [-1, 0.3 * q, 1]],
            energy_distribution_function='Single Level', energy_spread=0.1 * q)

# S_h1 = Trap('Shallow hole trap', charge_states=[[1, Eg - 0.1 * q, 1], [0, Eg - 0.1 * q, 1]],
#             energy_distribution_function='Single Level', energy_spread=0.01 * q)
# D_h1 = Trap('Deep hole trap', charge_states=[[1, Eg - 0.3 * q, 1], [0, Eg - 0.3 * q, 1]],
#             energy_distribution_function='Single Level', energy_spread=0.1 * q)

twist_dsl = Dislocation(b, 'Twist dislocation 1')
tilt_dsl = Dislocation(b, 'Tilt dislocation 1')
tilt_dsl.add_trap(S_e1, 9.0e7)  # Linear density of traps [1/m]
tilt_dsl.add_trap(D_e1, 5.0e7)
# tilt_dsl.add_trap(S_h1, 9.0e7) # Linear density of traps [1/m]
# tilt_dsl.add_trap(D_h1, 5.0e7)
print twist_dsl
print tilt_dsl

BW = BondingInterface(1.5e-7, 1.0e-7, 3.0, 0.5, twist_dsl, tilt_dsl)
# BW2 = BondingInterface(9.0e-7, 0.6e-7, 3.0, 0.5, twist_dsl, tilt_dsl)
print BW
Si.add_bonding_interface(BW)
# Si.add_bonding_interface(BW2)
'''

T_start = 75.0
T_stop = 80.0
T_step = 1.0

V_p = 3.0
V_rb = -2.0

T_range = np.linspace(T_start, T_stop, abs((T_start - T_stop) / T_step) + 1, endpoint=True)
points_number = T_range.size

Vd = np.zeros((T_range.size, 2))
Vd_err = np.zeros((T_range.size, 2))
J = np.zeros((T_range.size, 2))
J_err = np.zeros((T_range.size, 2))
Psi_list = [[] for i in range(T_range.size)]
Poisson_rel_err = 1e-1

def calculate(args):
    Si, Electrode, T, V_p = args
    sample_data_dir = join(dirname(__file__), '03-data')
    db_name = '03_Au_nSi_BW_transient' + '_%02.2fVp_%02.2fVrb_%03.2fK' % (V_p, abs(V_rb), T) + '.db'
    db_name = join(sample_data_dir, db_name)
    MyProject = Project(db_name=db_name, backend='sqlite', hostname='', overwrite=False)
    MyDiode = SchottkyDiode(MyProject, 'Au-Si_BW', Electrode, Si, DeadLayer=1.5e-7, L=5e-6)
    MyDiode.set_T(T)
    MyDiode.set_Va(V_p)
    print(T)

    type_sign = -1 if MyDiode.Semiconductor.dop_type == 'n' else 1
    Psi = Psi_approx(MyDiode.L, -(MyDiode.V_bi(eV=True) + type_sign * V_p), 0.0)

    Psi, E, z_nodes, rho_err_points, Vd, Vd_err, J, J_err, \
        BI_F, dopants_F, ic_id = Poisson.Reccurent_Poisson_solver(MyDiode, Psi, Vd_error=1e-6,
                                                                  equilibrium_filling=True, t=mp.inf,
                                                                  initial_condition_id=-1,
                                                                  rho_rel_err=Poisson_rel_err, max_iter=100,
                                                                  debug=False)
    return db_name


args = np.empty((len(T_range), 4), dtype=object)
args[:, 0] = silicon
args[:, 1] = electrode
args[:, 2] = T_range
args[:, 3] = V_p

pool = Pool()
db_names = np.array(pool.map(calculate, args))
print(db_names)

