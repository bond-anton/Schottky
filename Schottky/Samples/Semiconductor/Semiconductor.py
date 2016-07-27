from __future__ import division, print_function
import numbers

from Schottky.Samples import Sample
from Schottky.Samples.Semiconductor.Dopant import Dopant


class Semiconductor(Sample):

    def __init__(self, client, name, epsilon=None,
                 bands_density_of_states=None,
                 band_gap_parameters=None,
                 electron_mobility_parameters=None,
                 hole_mobility_parameters=None,
                 dopants=None, interfaces=None,
                 description=None):
        super(Semiconductor, self).__init__(client=client, name=name, description=description)
        self.load_create_sample()
        self.epsilon = None
        self.dopants = None
        self.interfaces = None
        self._read_in_epsilon(epsilon=epsilon)
        self._read_in_dopant(dopants)

    def _read_in_epsilon(self, epsilon):
        try:
            self.epsilon = self.parameters['epsilon'].float_value
        except KeyError:
            pass
        if self.epsilon != epsilon and epsilon is not None:
            self.set_epsilon(epsilon)

    def set_epsilon(self, epsilon):
        assert isinstance(epsilon, numbers.Number), 'epsilon must be a number'
        try:
            self.parameters['epsilon'].float_value = float(epsilon)
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_numeric_parameter(name='epsilon',
                                                                               value=float(epsilon),
                                                                               description='Permittivity of space')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.load_create_sample()
        self.epsilon = float(epsilon)

    def _read_in_dopant(self, dopant):
        try:
            trap_parameter = self.parameters['Trap']
            trap_module_name, trap_class_name, trap_name = trap_parameter.string_value.split('::')
            trap_id = int(trap_parameter.float_value)
            trap_module = __import__(trap_module_name, fromlist=[trap_class_name])
            trap_class = getattr(trap_module, trap_class_name)
            trap_sample = trap_class(client=self.client.session_manager, name=trap_name)
            if trap_sample.sample.id == trap_id:
                self.dopants = trap_sample
            else:
                print('Trap IDs do not match')
        except KeyError:
            pass
        if dopant is not None:
            if self.dopants != dopant:
                self.set_dopant(dopant)

    def set_dopant(self, dopant):
        assert isinstance(dopant, Dopant), 'Expected a Dopant instance'
        assert isinstance(dopant, Sample), 'Expected a Sample instance'
        string_value = dopant.__class__.__module__ + '::'
        string_value += dopant.__class__.__name__ + '::'
        string_value += dopant.name
        try:
            trap_parameter = self.parameters['Trap']
            trap_parameter.float_value = dopant.sample.id
            trap_parameter.string_value = string_value
            trap_parameter.description = dopant.description
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_generic_parameter(
                name='Dopant',
                float_value=float(dopant.sample.id),
                string_value=string_value,
                description=dopant.description)
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
        self.parameters = {}
        self.load_create_sample()
        self.dopants = dopant
