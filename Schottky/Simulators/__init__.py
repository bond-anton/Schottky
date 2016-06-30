from __future__ import division, print_function
import numpy as np

from ScientificProjects.Client import Client


class Simulator(object):

    def __init__(self, client, name=None, description=None, parts=None,
                 category_name=None, category_description=None):
        assert isinstance(client, Client), 'Valid ScientificProjects Client instance is required'
        self.client = client.user_manager
        if name is None:
            name = 'Software Simulator'
        self.name = str(name)
        if description is None:
            description = 'Simulates properties of corresponding entity'
        self.description = description
        self.parts = []
        if isinstance(parts, (list, tuple, np.ndarray)):
            for part in parts:
                if isinstance(part, Simulator):
                    self.parts.append(part)
        self.manufacturer = self._create_manufacturer()
        if category_name is None:
            category_name = 'Simulation'
        if category_description is None:
            category_description = 'Simulation software'
        self.category = self._create_categories(category_name, category_description)
        self.equipment = self._register_equipment()

    def _create_manufacturer(self):
        name = 'Anton Bondarenko'
        name_short = 'ABS'
        description = 'Scientific software tools'
        return self.client.equipment_manager.create_manufacturer(name=name, name_short=name_short,
                                                                 description=description)

    def _create_categories(self, category_name, category_description):
        description = 'Measurement, automation, control, simulation, and other software tools'
        root_category = self.client.equipment_manager.create_equipment_category(name='Software',
                                                                                description=description)
        return self.client.equipment_manager.create_equipment_category(name=category_name,
                                                                       description=category_description,
                                                                       parent=root_category)

    def _register_equipment(self):
        equipment = self.client.equipment_manager.create_equipment(name=self.name,
                                                                   category=self.category,
                                                                   manufacturer=self.manufacturer,
                                                                   description=self.description)
        if self.parts:
            assembly = self.client.equipment_manager.create_equipment_assembly(name=self.name+'-parts',
                                                                               description='Parts for ' + self.name)
            for part in self.parts:
                self.client.equipment_manager.add_equipment_to_assembly(assembly=assembly,
                                                                        equipment=part.equipment)
            equipment.assembly = assembly
        return equipment

    def __str__(self):
        description = 'Simulator: %s' % self.name
        description += '\n %s' % description
        description += '\n Category: %s' % self.category.name
        return description
