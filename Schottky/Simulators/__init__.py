from __future__ import division, print_function
import numpy as np

from ScientificProjects.Client import Client


class Simulator(object):

    def __init__(self, client, name=None, description=None, parts=None,
                 category=None, measurement_types=None, measurements=None):
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
        if category is None:
            category = {
                'name': 'Software',
                'description': 'Measurement, automation, control, simulation, and other software tools',
                'subcategory': {'name': 'Simulation',
                                'description': 'Simulation software',
                                'subcategory': None}}
        self.equipment = self._register_equipment(category)
        if measurement_types is None:
            measurement_types = [{'name': 'Synthetic data',
                                  'description': 'Pure synthetic data',
                                  'children': []}]
        if measurements is None:
            measurements = [{'name': 'random data generation',
                             'description': 'generates random data',
                             'type': 'Synthetic data'}]
        self.measurements = measurements
        self._register_measurement_types(measurement_types)

    def _create_manufacturer(self):
        name = 'Anton Bondarenko'
        name_short = 'ABS'
        description = 'Scientific software tools'
        return self.client.equipment_manager.create_manufacturer(name=name, name_short=name_short,
                                                                 description=description)

    def _create_categories(self, category, parent=None):
        if category is None:
            return parent
        parent = self.client.equipment_manager.create_equipment_category(name=category['name'],
                                                                         description=category['description'],
                                                                         parent=parent)
        return self._create_categories(category=category['subcategory'], parent=parent)

    def _register_equipment(self, category):
        category = self._create_categories(category=category)
        equipment = self.client.equipment_manager.create_equipment(name=self.name,
                                                                   category=category,
                                                                   manufacturer=self.manufacturer,
                                                                   description=self.description)
        if self.parts:
            assembly = self.client.equipment_manager.create_equipment_assembly(name=self.name + '-parts',
                                                                               description='Parts for ' + self.name)
            for part in self.parts:
                self.client.equipment_manager.add_equipment_to_assembly(assembly=assembly,
                                                                        equipment=part.equipment)
            equipment.assembly = assembly
        return equipment

    def _register_measurement_types(self, measurement_types, parent=None):
        for measurement_type in measurement_types:
            new_parent = self.client.measurement_type_manager.create_measurement_type(
                name=measurement_type['name'],
                description=measurement_type['description'],
                parent=parent)
            if measurement_type['children']:
                self._register_measurement_types(measurement_type['children'], parent=new_parent)
            else:
                self.client.equipment_manager.add_measurement_type_to_equipment(self.equipment, new_parent)

    def _register_measurement(self):
        print(self.name)

    def __str__(self):
        description = 'Simulator: %s' % self.name
        description += '\n %s' % description
        description += '\n Category: %s' % self.equipment.category.name
        return description
