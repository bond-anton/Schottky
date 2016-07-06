from __future__ import division
import numbers

from Schottky.Samples import Sample


class Trap(Sample):

    def __init__(self, client, name, description=None,
                 charge_state=None, activation_energy=None, band=None,
                 energy_distribution_function=None, energy_spread=None,
                 electron_capture_cross_section=None, electron_capture_cross_section_activation_energy=None,
                 hole_capture_cross_section=None, hole_capture_cross_section_activation_energy=None,
                 trap_potential=None):
        super(Trap, self).__init__(client=client, name=name, description=description)
        self.load_create_sample()
        self.charge_state = {}
        self.activation_energy = {}
        self.band = None
        self._read_in_charge_state(charge_state=charge_state)
        self._read_in_activation_energy(activation_energy=activation_energy)
        self._read_in_band(band=band)

    def _read_in_charge_state(self, charge_state):
        try:
            charge_states = self.parameters['Charge state'].children
            for parameter in charge_states:
                if parameter.name == 'empty':
                    self.charge_state['empty'] = int(parameter.float_value)
                elif parameter.name == 'full':
                    self.charge_state['full'] = int(parameter.float_value)
        except KeyError:
            pass
        if self.charge_state != charge_state:
            self.set_charge_state(charge_state)

    def set_charge_state(self, charge_state):
        assert isinstance(charge_state, dict), 'Charge state must be a dict with keys "empty" and "full"'
        assert 'empty' in charge_state.keys(), 'Charge state must be a dict with keys "empty" and "full"'
        assert 'full' in charge_state.keys(), 'Charge state must be a dict with keys "empty" and "full"'
        assert isinstance(charge_state['empty'], int), 'Charge state must be integer'
        assert isinstance(charge_state['full'], int), 'Charge state must be integer'
        try:
            charge_states = self.parameters['Charge state'].children
            for parameter in charge_states:
                if parameter.name == 'empty':
                    parameter.float_value = charge_state['empty']
                elif parameter.name == 'full':
                    parameter.float_value = charge_state['full']
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager._create_parameter(name='Charge state',
                                                                        description='Charge states of the trap')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='empty',
                                                                   value=charge_state['empty'],
                                                                   description='Charge state of empty trap',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='full',
                                                                   value=charge_state['full'],
                                                                   description='Charge state of full trap',
                                                                   parent=parameter)
            self.load_create_sample()
        self.charge_state = charge_state

    def _read_in_activation_energy(self, activation_energy):
        try:
            trap_energetics = self.parameters['Activation energy'].children
            for energy in trap_energetics:
                if energy.name == 'empty':
                    self.activation_energy['empty'] = energy.float_value
                elif energy.name == 'full':
                    self.activation_energy['full'] = energy.float_value
        except KeyError:
            pass
        if self.activation_energy != activation_energy:
            self.set_activation_energy(activation_energy)

    def set_activation_energy(self, activation_energy):
        assert isinstance(activation_energy, dict), 'Activation energy must be a dict with keys "empty" and "full"'
        assert 'empty' in activation_energy.keys(), 'Activation energy must be a dict with keys "empty" and "full"'
        assert 'full' in activation_energy.keys(), 'Activation energy must be a dict with keys "empty" and "full"'
        assert isinstance(activation_energy['empty'], numbers.Number), 'Activation energy must be a number'
        assert isinstance(activation_energy['full'], numbers.Number), 'Activation energy must be a number'
        try:
            energies = self.parameters['Activation energy'].children
            for parameter in energies:
                if parameter.name == 'empty':
                    parameter.float_value = float(activation_energy['empty'])
                elif parameter.name == 'full':
                    parameter.float_value = float(activation_energy['full'])
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager._create_parameter(name='Activation energy',
                                                                        description='Activation energy of the trap')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='empty',
                                                                   value=float(activation_energy['empty']),
                                                                   description='Activation energy of empty trap',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='full',
                                                                   value=float(activation_energy['full']),
                                                                   description='Activation energy of full trap',
                                                                   parent=parameter)
            self.load_create_sample()
        self.activation_energy = activation_energy

    def _read_in_band(self, band):
        try:
            self.band = self.parameters['Band'].string_value
        except KeyError:
            pass
        if self.band != band:
            self.set_band(band)

    def set_band(self, band):
        assert band in ['Ec', 'Ev'], 'Band must be either "Ec" or "Ev"'
        try:
            self.parameters['Band'].string_value = band
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_string_parameter(name='Band',
                                                                              value=band,
                                                                              description='Major band of the trap')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.load_create_sample()
        self.band = band
