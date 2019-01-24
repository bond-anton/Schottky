from __future__ import division, print_function

import numpy as np
from matplotlib import pyplot as plt

from BDMesh.Mesh1DUniform import Mesh1DUniform
from BDMesh.TreeMesh1DUniform import TreeMesh1DUniform

from Schottky.Dopant import Dopant
from Schottky.Semiconductor import Semiconductor
from Schottky.Reference import database
from Schottky import constant

reference = database[0]

c = Mesh1DUniform(0.0, 5e-6, physical_step=1e-6)
f = Mesh1DUniform(0.0, 5e-6, physical_step=1e-6)
c.solution = np.ones(c.num) * 1e21
f.solution = np.zeros(f.num)

silicon = Semiconductor('Si', reference)

phosphorus = Dopant('P', TreeMesh1DUniform(c, aligned=True), TreeMesh1DUniform(f, aligned=True),
                    0.045 * constant.q, silicon.band_gap_t(0.0) - 0.045 * constant.q,
                    1e-15, 1e-15)
phosphorus.charge_state = {0: +1, 1: 0}

silicon.dopants = [phosphorus]

t = np.linspace(0.0, 300.0, num=21, endpoint=True)
mu = np.linspace(0, silicon.band_gap_t(0.0), num=100)
fig, (ax1, ax2) = plt.subplots(2, 1)
ef = []
for t_i in t:
    ef.append(silicon.el_chem_pot(t_i) / constant.q)
    print('T:', t_i, 'K. Ef:', ef[-1], 'eV.')
    n_e = np.array([silicon.n_e_t(mu_i, t_i) for mu_i in mu])
    n_h = np.array([silicon.n_h_t(mu_i, t_i) for mu_i in mu])
    q_tot = np.array([silicon.bulk_charge(mu_i, t_i) for mu_i in mu])

    # ax.plot(mu, n_e, 'b-')
    # ax.plot(mu, n_h, 'r-')
    # ax.semilogy(mu / constant.q, abs(n_h - n_e), 'g-')
    ax1.semilogy(mu / constant.q, abs(q_tot), 'g-')
ax2.plot(t, np.array(ef), 'r-o')
plt.show()
