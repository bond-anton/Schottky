from __future__ import division, print_function
import numpy as np

from ScientificProjects.Client import Client


class Simulator(object):

    def __init__(self, client, parts=None):
        assert isinstance(client, Client), 'Valid ScientificProjects Client instance is required'
        self.client = client.user_manager
        self.name = 'Software Simulator'
        self.description = 'Simulates properties of corresponding entity'
        self.parts = []
        if isinstance(parts, (list, tuple, np.ndarray)):
            for part in parts:
                if isinstance(part, Simulator):
                    self.parts.append(part)
        category = 'Simulation'
        category_description = 'Simulation software'
        self.category = self.create_categories(category, category_description)

    def create_categories(self, category, category_description):
        description = 'Measurement, automation, control, simulation, and other software tools'
        root_category = self.client.equipment_manager.create_equipment_category(name='Software',
                                                                                description=description)
        return self.client.equipment_manager.create_equipment_category(name=category,
                                                                       description=category_description,
                                                                       parent=root_category)

    def register_equipment(self):
        print(self.name)

    def __str__(self):
        description = 'Simulator: %s' % self.name
        description += '\n %s' % description
        description += '\n Category: %s' % self.category.name
        return description
