from __future__ import division, print_function

import numpy as np

from BDMesh.Mesh1DUniform import Mesh1DUniform
from BDMesh.TreeMesh1DUniform import TreeMesh1DUniform

from Schottky.Dopant import Dopant
from Schottky.Semiconductor import Semiconductor
from Schottky.Reference import database
from Schottky import constant

reference = database[0]

c = Mesh1DUniform(0.0, 5e-6, physical_step=1e-6)
f = Mesh1DUniform(0.0, 5e-6, physical_step=1e-6)
c.solution = np.ones(c.num) * 1e15
f.solution = np.zeros(f.num)

silicon = Semiconductor('Si', reference)

phosphorus = Dopant('P', TreeMesh1DUniform(c, aligned=True), TreeMesh1DUniform(f, aligned=True),
                    0.045 * constant.q, silicon.band_gap_t(0.0) - 0.045 * constant.q,
                    1e-15, 1e-15)

silicon.dopants = [phosphorus]

temperature = 300.0
print(silicon.el_chem_pot(temperature) / constant.q)

