# author: anton

from __future__ import division
from os.path import dirname, join

import numpy as np

from Schottky.Notation import q
from Schottky.Metal import Metal
from Schottky.Semiconductor import Semiconductor, Trap, Dopant, Dislocation, BondingInterface
from Potential.Potential_3D import ConstantFieldPotential, SuperposedPotential, ChargedCylinderPotential, \
    DislocationDeformationPotential

data_dir = join(dirname(__file__), '03-data')
prefix = '03_Au_nSi_BW'

electrode = Metal('Au', q * 5.1)

dopant = Dopant('Phosphorus', 1.0e21, [[+1, 0.045 * q, 1], [0, 0.045 * q, 1]])

silicon = Semiconductor('Si', lookup=True)
silicon.add_dopant(dopant)

lattice_parameter = 5.43e-10
burgers_vector = lattice_parameter * np.sqrt(2) / 2

dp = DislocationDeformationPotential('Deformation', 5, 4e-10)
cc = ChargedCylinderPotential('Charged Dislocation', charge_sign=-1, linear_charge_density=1e7 * 0,
                              radius=1e-9, epsilon=11.8)
ef = ConstantFieldPotential('External Field', (0.0, 0.0, 0.0))
sp = SuperposedPotential('Superposed', [dp, cc, ef])

shallow_electron_trap = Trap('Shallow electron trap', charge_states=[[0, 0.15 * q, 1], [-1, 0.15 * q, 1]],
                             energy_distribution_function='Single Level', energy_spread=0.01 * q, trap_potential=sp)
deep_electron_trap = Trap('Deep electron trap', charge_states=[[0, 0.3 * q, 1], [-1, 0.3 * q, 1]],
                          energy_distribution_function='Single Level', energy_spread=0.1 * q, trap_potential=sp)

twist_dislocation = Dislocation(burgers_vector, 'Twist dislocation 1')
tilt_dislocation = Dislocation(burgers_vector, 'Tilt dislocation 1')
tilt_dislocation.add_trap(shallow_electron_trap, 5.0e7)  # Linear density of traps [1/m]
tilt_dislocation.add_trap(deep_electron_trap, 1.0e7)

print twist_dislocation
print tilt_dislocation

bonding_interface = BondingInterface(1.5e-7, 1.0e-7, 3.0, 0.5, twist_dislocation, tilt_dislocation)
print bonding_interface

silicon.add_bonding_interface(bonding_interface)
