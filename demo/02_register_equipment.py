from __future__ import division, print_function

from ScientificProjects.Client import Client

from Schottky.Metal import Metal
from Schottky.Semiconductor import Trap, Dopant, Semiconductor
from Schottky import Diode
from Schottky._version import __version__

client = Client(config_file_name='config.ini')

client.user_manager.sign_in('bond_anton', 'secret_password')

description = 'My home brew software for measurements and simulations'
manufacturer = client.user_manager.equipment_manager.create_manufacturer(name='Anton Bondarenko Software',
                                                                         name_short='ABS',
                                                                         description=description)

description = 'Measurement and simulation software'
category = client.user_manager.equipment_manager.create_equipment_category(name='Software',
                                                                           description=description)
description = 'Simulation software'
category = client.user_manager.equipment_manager.create_equipment_category(name='Simulation Software',
                                                                           description=description,
                                                                           parent=category)

metal_simulator = client.user_manager.equipment_manager.create_equipment(name=Metal.equipment,
                                                                         description=Metal.description,
                                                                         category=category,
                                                                         manufacturer=manufacturer)
trap_simulator = client.user_manager.equipment_manager.create_equipment(name=Trap.equipment,
                                                                        description=Trap.description,
                                                                        category=category,
                                                                        manufacturer=manufacturer)

dopant_simulator = client.user_manager.equipment_manager.create_equipment(name=Dopant.equipment,
                                                                          description=Dopant.description,
                                                                          category=category,
                                                                          manufacturer=manufacturer)

ssa_name = 'Semiconductor simulator tools'
ssa_description = 'Semiconductor simulator incorporates many different tools inside'
semiconductor_assembly = client.user_manager.equipment_manager.create_equipment_assembly(name=ssa_name,
                                                                                         description=ssa_description)
client.user_manager.equipment_manager.add_equipment_to_assembly(semiconductor_assembly, trap_simulator)
client.user_manager.equipment_manager.add_equipment_to_assembly(semiconductor_assembly, dopant_simulator)
semiconductor_simulator = client.user_manager.equipment_manager.create_equipment(name=Semiconductor.equipment,
                                                                                 description=Semiconductor.description,
                                                                                 assembly=semiconductor_assembly,
                                                                                 category=category,
                                                                                 manufacturer=manufacturer)
print(semiconductor_simulator)

client.user_manager.sign_out()
