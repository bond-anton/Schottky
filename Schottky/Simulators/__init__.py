from __future__ import division, print_function
import numpy as np

from BDProjects.Client import Client

from Schottky.Samples import Sample


class Simulator(object):

    def __init__(self, client, name=None, description=None, label=None,
                 samples=None, parts=None,
                 category=None, measurement_types=None, measurements=None):
        assert isinstance(client, Client), 'Valid ScientificProjects Client instance is required.'
        assert client.user_manager.project_manager.project_opened(), 'Please open project to work with.'
        self.client = client.user_manager
        if name is None:
            name = 'Software Simulator'
        self.name = str(name)
        if label is None:
            label = self.name
        self.label = str(label)
        if description is None:
            description = 'Simulates properties of corresponding entity'
        self.description = description
        self.samples = {}
        if isinstance(samples, (list, tuple, np.ndarray)):
            for sample in samples:
                if isinstance(sample, Sample):
                    self.samples[sample.name] = sample
        self.parts = {}
        if isinstance(parts, (list, tuple, np.ndarray)):
            for part in parts:
                if isinstance(part, Simulator):
                    self.parts[part.label] = part
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
            measurement_types = [{'name': 'Simulated data',
                                  'description': 'Software simulated data',
                                  'children': []}]
        if measurements is None:
            measurements = [{'name': 'random data generation',
                             'description': 'generates random data',
                             'type': 'Simulated data'}]
        self.measurements = measurements
        self._register_measurement_types(measurement_types)

    def _create_manufacturer(self):
        name = 'BondDevices'
        name_short = 'BD'
        description = 'BondDevices Scientific Services'
        manufacturers = self.client.equipment_manager.get_manufacturers(name=name)
        if len(manufacturers) == 1:
            return manufacturers[0]
        elif len(manufacturers) == 0:
            return self.client.equipment_manager.create_manufacturer(name=name, name_short=name_short,
                                                                     description=description)
        else:
            raise ValueError('More than one manufacturer found for given name. Check the database.')

    def _create_categories(self, category, parent=None):
        if category is None:
            return parent
        categories = self.client.equipment_manager.get_equipment_categories(name=category['name'])
        if len(categories) == 1 and categories[0].parent == parent:
            parent = categories[0]
        elif len(categories) == 0:
            parent = self.client.equipment_manager.create_equipment_category(name=category['name'],
                                                                             description=category['description'],
                                                                             parent=parent)
        else:
            found = False
            for item in categories:
                if item.parent == parent:
                    parent = item
                    found = True
                    break
            if not found:
                parent = self.client.equipment_manager.create_equipment_category(name=category['name'],
                                                                                 description=category['description'],
                                                                                 parent=parent)
        return self._create_categories(category=category['subcategory'], parent=parent)

    def _register_equipment(self, category):
        category = self._create_categories(category=category)
        equipment_found = self.client.equipment_manager.get_equipment(name=self.name)
        if len(equipment_found) == 1 and equipment_found[0].category.name == category.name:
            equipment = equipment_found[0]
        elif len(equipment_found) == 0:
            equipment = self.client.equipment_manager.create_equipment(name=self.name,
                                                                       category=category,
                                                                       manufacturer=self.manufacturer,
                                                                       description=self.description)
        else:
            raise ValueError('More than one equipment items found for given name')
        if self.parts:
            assemblies = self.client.equipment_manager.get_equipment_assembly(name=self.name + '-parts')
            if len(assemblies) == 1:
                assembly = assemblies[0]
            elif len(assemblies) == 0:
                assembly = self.client.equipment_manager.create_equipment_assembly(name=self.name + '-parts',
                                                                                   description='Parts for ' + self.name)
            else:
                raise ValueError('More than one assembly found. Check database.')
            for part in self.parts.values():
                if part.equipment not in assembly.parts:
                    self.client.equipment_manager.add_equipment_to_assembly(assembly=assembly,
                                                                            equipment=part.equipment)
            equipment.assembly = assembly
        return equipment

    def _register_measurement_types(self, measurement_types, parent=None):
        for measurement_type in measurement_types:
            found = False
            types_found = self.client.measurement_type_manager.get_measurement_types(name=measurement_type['name'])
            if len(types_found) == 1:
                new_parent = types_found[0]
                found = True
            elif len(types_found) == 0:
                new_parent = self.client.measurement_type_manager.create_measurement_type(
                    name=measurement_type['name'],
                    description=measurement_type['description'],
                    parent=parent)
            else:
                for item in types_found:
                    if item.name == measurement_type['name']:
                        found = True
                        new_parent = item
                        break
                if not found:
                    new_parent = self.client.measurement_type_manager.create_measurement_type(
                        name=measurement_type['name'],
                        description=measurement_type['description'],
                        parent=parent)
            if measurement_type['children']:
                self._register_measurement_types(measurement_type['children'], parent=new_parent)
            else:
                if not found:
                    self.client.equipment_manager.add_measurement_type_to_equipment(self.equipment, new_parent)

    def measurement_lookup(self, measurement_details, parameters=None, input_data=None, progress_threshold=0):
        if parameters is None:
            parameters = []
        samples_parameters = []
        for sample in self.samples.values():
            samples_parameters += sample.parameters.values()
        measurement_parameters = parameters + samples_parameters
        measurements = self.client.measurement_manager.get_measurements(name=measurement_details['name'])
        matched_measurements = []
        for measurement in measurements:
            # check progress threshold
            if measurement.progress < progress_threshold:
                continue
            # check input_data match
            if measurement.input_data != input_data:
                continue
            # check measurement parameters match
            if len(measurement.parameters) == len(measurement_parameters):
                test_parameters = [] + measurement_parameters
                for parameter in measurement.parameters:
                    parameters_match = False
                    for i in range(len(test_parameters)):
                        if parameter.equals(test_parameters[i]):
                            parameters_match = True
                            del test_parameters[i]
                            break
                    if not parameters_match:
                        break
                if not parameters_match:
                    continue
            # check samples match
            samples = [i.sample for i in self.samples.values()]
            if len(measurement.samples) == len(samples):
                samples_match = True
                for sample in measurement.samples:
                    if sample not in samples:
                        samples_match = False
                        break
                if not samples_match:
                    continue
            matched_measurements.append(measurement)
        return matched_measurements

    def measurement_new(self, measurement_details, parameters=None, input_data=None):
        if parameters is None:
            parameters = []
        samples_parameters = []
        for sample in self.samples.values():
            samples_parameters += sample.parameters.values()
        measurement_parameters = parameters + samples_parameters
        measurement = self.client.measurement_manager.create_measurement(
            name=measurement_details['name'],
            measurement_type=measurement_details['type'],
            equipment=self.equipment,
            description=measurement_details['description'])
        if samples_parameters:
            copied_parameters = [self.client.parameter_manager.copy_parameter(p) for p in samples_parameters]
            for parameter in copied_parameters:
                self.client.measurement_manager.add_parameter_to_measurement(
                    measurement=measurement,
                    parameter=parameter)
        if parameters:
            for parameter in parameters:
                self.client.measurement_manager.add_parameter_to_measurement(
                    measurement=measurement,
                    parameter=parameter)
        for sample in self.samples.values():
            self.client.measurement_manager.add_sample_to_measurement(
                measurement=measurement,
                sample=sample.sample)
        if input_data is not None:
            self.client.measurement_manager.add_input_data_to_measurement(
                measurement=measurement,
                measurements_collection=input_data)
        return measurement

    def load_create_data_channel(self, channel_name, measurement, description, unit_name=None):
        channels = self.client.measurement_manager.get_data_channels(measurement=measurement, name=channel_name,
                                                                     exact=True)
        if len(channels) == 1:
            return channels[0]
        elif len(channels) == 0:
            return self.client.measurement_manager.create_data_channel(name=channel_name, measurement=measurement,
                                                                       description=description, unit_name=unit_name)
        else:
            raise ValueError('More than one channel found for given name (%s). Check the database.' % channel_name)

    def __str__(self):
        description = 'Simulator: %s' % self.name
        description += '\n %s' % description
        description += '\n Category: %s' % self.equipment.category.name
        return description
