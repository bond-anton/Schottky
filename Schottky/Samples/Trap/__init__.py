from __future__ import division
import numbers

from Schottky.Samples import Sample

energy_distribution_functions = {'Single Level': ['single level', 'monoenergetic level', 'single', 'monoenergetic'],
                                 'Gaussian Level': ['gaussian level', 'gaussian'],
                                 'Rectangular Level': ['rectangular level', 'rectangular', 'box']}


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
        self.energy_distribution_function = None
        self.energy_spread = None
        self.electron_capture_cross_section = None
        self.electron_capture_cross_section_activation_energy = None
        self.hole_capture_cross_section = None
        self.hole_capture_cross_section_activation_energy = None
        self.trap_potential = None
        self._read_in_charge_state(charge_state=charge_state)
        self._read_in_activation_energy(activation_energy=activation_energy)
        self._read_in_band(band=band)
        self._read_in_energy_distribution_function(energy_distribution_function=energy_distribution_function)
        self._read_in_energy_spread(energy_spread=energy_spread)
        self._read_in_electron_capture_cross_section(capture_cross_section=electron_capture_cross_section)
        self._read_in_electron_capture_cross_section_activation_energy(electron_capture_cross_section_activation_energy)
        self._read_in_hole_capture_cross_section(capture_cross_section=hole_capture_cross_section)
        self._read_in_hole_capture_cross_section_activation_energy(hole_capture_cross_section_activation_energy)
        self.trap_potential = trap_potential

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
        if self.charge_state != charge_state and charge_state is not None:
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
        if self.activation_energy != activation_energy and activation_energy is not None:
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
        if self.band != band and band is not None:
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

    def _read_in_energy_distribution_function(self, energy_distribution_function):
        try:
            self.energy_distribution_function = self.parameters['Energy distribution function'].string_value
        except KeyError:
            pass
        if energy_distribution_function is not None:
            for key in energy_distribution_functions.keys():
                if energy_distribution_function.lower() in energy_distribution_functions[key]:
                    energy_distribution_function = key
                    break
            if self.energy_distribution_function != energy_distribution_function:
                self.set_energy_distribution_function(energy_distribution_function)

    def set_energy_distribution_function(self, energy_distribution_function):
        variants = [j for i in energy_distribution_functions.keys() for j in energy_distribution_functions[i]]
        assert energy_distribution_function in variants, 'Energy distribution function must be a defined string'
        for key in energy_distribution_functions.keys():
            if energy_distribution_function.lower() in energy_distribution_functions[key]:
                energy_distribution_function = key
                break
        try:
            self.parameters['Energy distribution function'].string_value = energy_distribution_function
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_string_parameter(
                name='Energy distribution function',
                value=energy_distribution_function,
                description='Energy distribution function of the trap')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.load_create_sample()
        self.band = energy_distribution_function

    def _read_in_energy_spread(self, energy_spread):
        try:
            self.energy_spread = self.parameters['Energy spread'].float_value
        except KeyError:
            pass
        if self.energy_spread != energy_spread and energy_spread is not None:
            self.set_energy_spread(energy_spread)

    def set_energy_spread(self, energy_spread):
        assert isinstance(energy_spread, numbers.Number), 'Energy spread must be a number'
        try:
            self.parameters['Energy spread'].float_value = energy_spread
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_numeric_parameter(name='Energy spread',
                                                                               value=energy_spread,
                                                                               unit_name='eV',
                                                                               description='Energy spread of the trap')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.load_create_sample()
        self.energy_spread = energy_spread

    def _read_in_electron_capture_cross_section(self, capture_cross_section):
        try:
            self.electron_capture_cross_section = self.parameters['Electron capture cross section'].float_value
        except KeyError:
            pass
        if self.electron_capture_cross_section != capture_cross_section and capture_cross_section is not None:
            self.set_electron_capture_cross_section(capture_cross_section)

    def set_electron_capture_cross_section(self, capture_cross_section):
        assert isinstance(capture_cross_section, numbers.Number), 'Capture cross section must be a number'
        try:
            self.parameters['Electron capture cross section'].float_value = capture_cross_section
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_numeric_parameter(
                name='Electron capture cross section',
                value=capture_cross_section,
                unit_name='cm^-2',
                description='Electron capture cross section of the trap')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.load_create_sample()
        self.electron_capture_cross_section = capture_cross_section

    def _read_in_electron_capture_cross_section_activation_energy(self, activation_energy):
        try:
            parameter_value = self.parameters['Electron capture cross section activation energy'].float_value
            self.electron_capture_cross_section_activation_energy = parameter_value
        except KeyError:
            pass
        if self.electron_capture_cross_section_activation_energy != activation_energy and activation_energy is not None:
            self.set_electron_capture_cross_section_activation_energy(activation_energy)

    def set_electron_capture_cross_section_activation_energy(self, activation_energy):
        assert isinstance(activation_energy, numbers.Number), 'Activation energy must be a number'
        try:
            self.parameters['Electron capture cross section activation energy'].float_value = activation_energy
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_numeric_parameter(
                name='Electron capture cross section activation energy',
                value=activation_energy,
                unit_name='eV',
                description='Electron capture cross section activation energy of the trap')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.load_create_sample()
        self.electron_capture_cross_section_activation_energy = activation_energy

    def _read_in_hole_capture_cross_section(self, capture_cross_section):
        try:
            self.hole_capture_cross_section = self.parameters['Hole capture cross section'].float_value
        except KeyError:
            pass
        if self.hole_capture_cross_section != capture_cross_section and capture_cross_section is not None:
            self.set_hole_capture_cross_section(capture_cross_section)

    def set_hole_capture_cross_section(self, capture_cross_section):
        assert isinstance(capture_cross_section, numbers.Number), 'Capture cross section must be a number'
        try:
            self.parameters['Hole capture cross section'].float_value = capture_cross_section
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_numeric_parameter(
                name='Hole capture cross section',
                value=capture_cross_section,
                unit_name='cm^-2',
                description='Hole capture cross section of the trap')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.load_create_sample()
        self.hole_capture_cross_section = capture_cross_section

    def _read_in_hole_capture_cross_section_activation_energy(self, activation_energy):
        try:
            parameter_value = self.parameters['Hole capture cross section activation energy'].float_value
            self.hole_capture_cross_section_activation_energy = parameter_value
        except KeyError:
            pass
        if self.hole_capture_cross_section_activation_energy != activation_energy and activation_energy is not None:
            self.set_hole_capture_cross_section_activation_energy(activation_energy)

    def set_hole_capture_cross_section_activation_energy(self, activation_energy):
        assert isinstance(activation_energy, numbers.Number), 'Activation energy must be a number'
        try:
            self.parameters['Hole capture cross section activation energy'].float_value = activation_energy
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_numeric_parameter(
                name='Hole capture cross section activation energy',
                value=activation_energy,
                unit_name='eV',
                description='Hole capture cross section activation energy of the trap')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.load_create_sample()
        self.hole_capture_cross_section_activation_energy = activation_energy