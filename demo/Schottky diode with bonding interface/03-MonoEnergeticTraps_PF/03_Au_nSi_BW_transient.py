# -*- coding: utf-8 -*-

'''
Created on 10 апр. 2015 г.

@author: anton
'''

from __future__ import division
from os.path import dirname, join

from pathos.pools import ProcessPool as Pool

import numpy as np
import pandas as pd
np.seterr(all='raise')
import mpmath as mp
# from scipy.interpolate import interp1d
# mp.mp.dps = 25

#import matplotlib
# matplotlib.use("Qt4Agg")
#from matplotlib import pyplot as plt

# plt.xkcd()

from diode_description import *

from ProjectManager import Project

#from Schottky.Notation import q, Eg
#from Schottky.Metal import Metal
#from Schottky.Semiconductor import Semiconductor, Trap, Dopant, Dislocation, BondingInterface
from Schottky.Diode import SchottkyDiode, Poisson, Kinetics, Visual
from Schottky.Helpers import Psi_approx

colors = ['b', 'g', 'y', 'k', 'm', 'c', 'b', 'g', 'y', 'k', 'm', 'c', 'b', 'g', 'y', 'k', 'm', 'c', 'b', 'g', 'y', 'k',
          'm', 'c']

T_start = 41.0
T_stop = 41.0
T_step = 1.0

V_p = 1.0
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
    Si, Electrode, T, V_p, V_rb = args
    sample_data_dir = join(dirname(__file__), '03-data')
    db_name = '03_Au_nSi_BW_transient' + '_%02.2fVp_%02.2fVrb_%03.2fK' % (V_p, abs(V_rb), T) + '.db'
    db_name = join(sample_data_dir, db_name)
    MyProject = Project(db_name=db_name, backend='sqlite', hostname='', overwrite=False)
    MyDiode = SchottkyDiode(MyProject, 'Au-Si_BW', Electrode, Si, DeadLayer=1.5e-7, L=5e-6)
    MyDiode.set_T(T)
    MyDiode.set_Va(V_p)

    type_sign = -1 if MyDiode.Semiconductor.dop_type == 'n' else 1
    Psi = Psi_approx(MyDiode.L, -(MyDiode.V_bi(eV=True) + type_sign * V_p), 0.0)

    Psi, E, z_nodes, rho_err_points, Vd, Vd_err, J, J_err, \
        BI_F, dopants_F, ic_id = Poisson.Reccurent_Poisson_solver(MyDiode, Psi, Vd_error=1e-6,
                                                                  equilibrium_filling=True, t=mp.inf,
                                                                  initial_condition_id=-1,
                                                                  rho_rel_err=Poisson_rel_err, max_iter=100,
                                                                  debug=False)

    #plt.plot(z_nodes, -Psi(z_nodes), 'r-o')
    #plt.show()

    MyDiode.set_Va(V_rb)
    Psi, E, z_nodes, rho_err_points, Vd, Vd_err, J, J_err, \
        BI_F, dopants_F, ic_id = Poisson.Reccurent_Poisson_solver(MyDiode, Psi, Vd_error=1e-6,
                                                                  equilibrium_filling=False, fast_traps=['Phosphorus'],
                                                                  t=0.0,
                                                                  initial_condition_id=ic_id,
                                                                  rho_rel_err=Poisson_rel_err, max_iter=100,
                                                                  debug=False)

    #plt.plot(z_nodes, -Psi(z_nodes), 'r-o')
    #plt.show()

    t_points, potential_t, field_d, z_t, diode_voltage_drop_t, current_density_t, \
        bonding_interfaces_f_t, dopants_f_t, last_state_id = Kinetics.traps_kinetics(MyDiode, ic_id,
                                                                                     1e-9, 10e-3, 100e-3,
                                                                                     fast_traps=['Phosphorus'],
                                                                                     rho_rel_err=Poisson_rel_err,
                                                                                     df_threshold=1e-2, debug=True)

    return [T, t_points, bonding_interfaces_f_t, dopants_f_t]


args = np.empty((len(T_range), 5), dtype=object)
args[:, 0] = silicon
args[:, 1] = electrode
args[:, 2] = T_range
args[:, 3] = V_p
args[:, 4] = V_rb

pool = Pool()
results = np.array(pool.map(calculate, args))
print results

sample_data_dir = join(dirname(__file__), '03-kinetics')

t1 = np.array([1e-6, 1e-5, 2e-5, 5e-5, 1e-5])
dlts_array = np.zeros((len(T_range), len(t1)), dtype=np.float)
#_, ax = plt.subplots()
for result in results:
    T, t_points, bonding_interfaces_f_t, dopants_f_t = result
    csv_name = '03_Au_nSi_BW_transient' + '_%02.2fVp_%02.2fVrb_%03.2fK' % (V_p, abs(V_rb), T) + '.csv'
    csv_name = join(sample_data_dir, csv_name)
    df = pd.DataFrame(bonding_interfaces_f_t)
    df['time'] = t_points
    df.set_index('time', inplace=True)
    df.to_csv(csv_name)
#ax.legend_.remove()
#plt.show()