from __future__ import division, print_function
import numbers

from Schottky.Samples import Sample
from Schottky.Samples.Trap import Trap


class Dopant(Sample):

    def __init__(self, client, name, trap=None, concentration=None, description=None):
        super(Dopant, self).__init__(client=client, name=name, description=description)
        self.load_create_sample()
        self.concentration = None
        self.trap = None
        self._read_in_concentration(concentration)
        self._read_in_trap(trap)

    def _read_in_concentration(self, concentration):
        try:
            self.concentration = self.parameters['Concentration'].float_value
        except KeyError:
            pass
        if self.concentration != concentration and concentration is not None:
            self.set_concentration(concentration)

    def set_concentration(self, concentration):
        assert isinstance(concentration, numbers.Number), 'Concentration must be a number'
        try:
            self.parameters['Concentration'].float_value = concentration
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_numeric_parameter(name='Concentration',
                                                                               value=concentration,
                                                                               unit_name='m^-3',
                                                                               description='Dopant net concentration')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.load_create_sample()
        self.concentration = concentration

    def _read_in_trap(self, trap):
        try:
            trap_parameter = self.parameters['Trap']
            trap_module_name, trap_class_name, trap_name = trap_parameter.string_value.split('::')
            trap_id = int(trap_parameter.float_value)
            trap_module = __import__(trap_module_name, fromlist=[str(trap_class_name)])
            trap_class = getattr(trap_module, trap_class_name)
            trap_sample = trap_class(client=self.client.session_manager, name=trap_name)
            if trap_sample.sample.id == trap_id:
                self.trap = trap_sample
            else:
                print('Trap IDs do not match')
        except KeyError:
            pass
        if trap is not None:
            if self.trap != trap:
                self.set_trap(trap)

    def set_trap(self, trap):
        assert isinstance(trap, Trap), 'Expected a Trap instance'
        assert isinstance(trap, Sample), 'Expected a Sample instance'
        string_value = trap.__class__.__module__ + '::'
        string_value += trap.__class__.__name__ + '::'
        string_value += trap.name
        try:
            trap_parameter = self.parameters['Trap']
            trap_parameter.float_value = trap.sample.id
            trap_parameter.string_value = string_value
            trap_parameter.description = trap.description
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_generic_parameter(
                name='Trap',
                float_value=float(trap.sample.id),
                string_value=string_value,
                description=trap.description)
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
        self.parameters = {}
        self.load_create_sample()
        self.trap = trap
