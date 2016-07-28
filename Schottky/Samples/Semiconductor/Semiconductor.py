from __future__ import division, print_function
import numbers

from Schottky.Samples import Sample
from Schottky.Samples.Semiconductor.Dopant import Dopant
from Schottky.Samples.Semiconductor.Interface import Interface


class Semiconductor(Sample):

    def __init__(self, client, name,
                 epsilon=None,
                 affinity=None,
                 effective_mass=None,
                 bands_density_of_states=None,
                 band_gap_parameters=None,
                 electron_mobility_parameters=None,
                 hole_mobility_parameters=None,
                 thermo_emission_coefficient=None,
                 dopants=None, interfaces=None,
                 description=None):
        super(Semiconductor, self).__init__(client=client, name=name, description=description)
        self.load_create_sample()
        self.epsilon = None
        self.affinity = None
        self.effective_mass = {}
        self.bands_density_of_states = {}
        self.band_gap_parameters = {}
        self.electron_mobility_parameters = {}
        self.hole_mobility_parameters = {}
        self.thermo_emission_coefficient = {}
        self.dopants = []
        self.interfaces = []
        self._read_in_epsilon(epsilon)
        self._read_in_affinity(affinity)
        self._read_in_effective_mass(effective_mass)
        self._read_in_bands_density_of_states(bands_density_of_states)
        self._read_in_band_gap_parameters(band_gap_parameters)
        self._read_in_electron_mobility_parameters(electron_mobility_parameters)
        self._read_in_hole_mobility_parameters(hole_mobility_parameters)
        self._read_in_thermo_emission_coefficient(thermo_emission_coefficient)
        self._read_in_dopants(dopants)
        self._read_in_interfaces(interfaces)

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

    def _read_in_affinity(self, affinity):
        try:
            self.affinity = self.parameters['affinity'].float_value
        except KeyError:
            pass
        if self.affinity != affinity and affinity is not None:
            self.set_affinity(affinity)

    def set_affinity(self, affinity):
        assert isinstance(affinity, numbers.Number), 'affinity must be a number'
        try:
            self.parameters['affinity'].float_value = float(affinity)
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_numeric_parameter(name='affinity',
                                                                               value=float(affinity),
                                                                               description='Electron affinity')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.load_create_sample()
        self.affinity = float(affinity)

    def _read_in_effective_mass(self, effective_mass):
        try:
            effective_mass_parameters = self.parameters['Effective mass'].children
            for particle in effective_mass_parameters:
                if particle.name == 'electron':
                    self.effective_mass['electron'] = particle.float_value
                elif particle.name == 'hole':
                    self.effective_mass['hole'] = particle.float_value
        except KeyError:
            pass
        if effective_mass is not None:
            if self.effective_mass != effective_mass:
                self.set_effective_mass(effective_mass)

    def set_effective_mass(self, effective_mass):
        assert isinstance(effective_mass, dict), 'Effective mass must be dictionary'
        assert 'electron' in effective_mass.keys(), 'Effective mass must be dictionary with electron and hole'
        assert 'hole' in effective_mass.keys(), 'Effective mass must be dictionary with electron and hole'
        assert isinstance(effective_mass['electron'], numbers.Number), 'Effective mass must be number'
        assert isinstance(effective_mass['hole'], numbers.Number), 'Effective mass must be number'
        try:
            effective_mass_parameters = self.parameters['Effective mass'].children
            for particle in effective_mass_parameters:
                if particle.name == 'electron':
                    particle.float_value = float(effective_mass['electron'])
                elif particle.name == 'hole':
                    particle.float_value = float(effective_mass['hole'])
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_dict_parameter(name='Effective mass',
                                                                            description='Effective masses of carriers')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='electron',
                                                                   value=float(effective_mass['electron']),
                                                                   description='Electron effective mass',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='hole',
                                                                   value=float(effective_mass['hole']),
                                                                   description='Hole effective mass',
                                                                   parent=parameter)
            self.load_create_sample()

    def _read_in_bands_density_of_states(self, bands_density_of_states):
        try:
            bands_parameters = self.parameters['Bands density of states'].children
            for band in bands_parameters:
                if band.name == 'Nc':
                    self.bands_density_of_states['Nc'] = band.float_value
                elif band.name == 'Nv':
                    self.bands_density_of_states['Nv'] = band.float_value
        except KeyError:
            pass
        if bands_density_of_states is not None:
            if self.bands_density_of_states != bands_density_of_states:
                self.set_bands_density_of_states(bands_density_of_states)

    def set_bands_density_of_states(self, bands_density_of_states):
        assert isinstance(bands_density_of_states, dict), 'Bands density of states must be dictionary'
        assert 'Nc' in bands_density_of_states.keys(), 'Bands density of states must be dictionary with Nc and Nv value'
        assert 'Nv' in bands_density_of_states.keys(), 'Bands density of states must be dictionary with Nc and Nv value'
        assert isinstance(bands_density_of_states['Nc'], numbers.Number), 'Band density of states must be number'
        assert isinstance(bands_density_of_states['Nv'], numbers.Number), 'Band density of states must be number'
        try:
            bands_parameters = self.parameters['Bands density of states'].children
            for band in bands_parameters:
                if band.name == 'Nc':
                    band.float_value = float(bands_density_of_states['Nc'])
                elif band.name == 'Nv':
                    band.float_value = float(bands_density_of_states['Nv'])
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_dict_parameter(name='Bands density of states',
                                                                            description='Bands density of states')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='Nc',
                                                                   value=float(bands_density_of_states['Nc']),
                                                                   description='Nc density of states',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='Nv',
                                                                   value=float(bands_density_of_states['Nv']),
                                                                   description='Nv density of states',
                                                                   parent=parameter)
            self.load_create_sample()

    def _read_in_band_gap_parameters(self, band_gap_parameters):
        try:
            band_gap_parameter = self.parameters['Band gap parameters'].children
            for component in band_gap_parameter:
                if component.name == 'Eg_0':
                    self.band_gap_parameters['Eg_0'] = component.float_value
                elif component.name == 'alpha':
                    self.band_gap_parameters['alpha'] = component.float_value
                elif component.name == 'beta':
                    self.band_gap_parameters['beta'] = component.float_value
        except KeyError:
            pass
        if band_gap_parameters is not None:
            if self.band_gap_parameters != band_gap_parameters:
                self.set_band_gap_parameters(band_gap_parameters)

    def set_band_gap_parameters(self, band_gap_parameters):
        assert isinstance(band_gap_parameters, dict), 'Band gap parameters must be dictionary'
        assert 'Eg_0' in band_gap_parameters.keys(), 'Band gap parameters must contain Eg_0, alpha, and beta'
        assert 'alpha' in band_gap_parameters.keys(), 'Band gap parameters must contain Eg_0, alpha, and beta'
        assert 'beta' in band_gap_parameters.keys(), 'Band gap parameters must contain Eg_0, alpha, and beta'
        assert isinstance(band_gap_parameters['Eg_0'], numbers.Number), 'Eg_0 must be number'
        assert isinstance(band_gap_parameters['alpha'], numbers.Number), 'alpha must be number'
        assert isinstance(band_gap_parameters['beta'], numbers.Number), 'beta must be number'
        try:
            band_gap_parameter = self.parameters['Band gap parameters'].children
            for component in band_gap_parameter:
                if component.name == 'Eg_0':
                    component.float_value = float(band_gap_parameters['Eg_0'])
                elif component.name == 'alpha':
                    component.float_value = float(band_gap_parameters['alpha'])
                elif component.name == 'beta':
                    component.float_value = float(band_gap_parameters['beta'])
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_dict_parameter(name='Band gap parameters',
                                                                            description='Band gap parameters')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='Eg_0',
                                                                   value=float(band_gap_parameters['Eg_0']),
                                                                   description='Eg_0',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='alpha',
                                                                   value=float(band_gap_parameters['alpha']),
                                                                   description='alpha parameter',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='beta',
                                                                   value=float(band_gap_parameters['beta']),
                                                                   description='beta parameter',
                                                                   parent=parameter)
            self.load_create_sample()

    def _read_in_electron_mobility_parameters(self, electron_mobility_parameters):
        try:
            mobility_parameter = self.parameters['Electron mobility parameters'].children
            for component in mobility_parameter:
                if component.name == 'mu_L0':
                    self.electron_mobility_parameters['mu_L0'] = component.float_value
                elif component.name == 'v_s':
                    self.electron_mobility_parameters['v_s'] = component.float_value
                elif component.name == 'A':
                    self.electron_mobility_parameters['A'] = component.float_value
                elif component.name == 'B':
                    self.electron_mobility_parameters['B'] = component.float_value
                elif component.name == 'alpha':
                    self.electron_mobility_parameters['alpha'] = component.float_value
                elif component.name == 'beta':
                    self.electron_mobility_parameters['beta'] = component.float_value
        except KeyError:
            pass
        if electron_mobility_parameters is not None:
            if self.electron_mobility_parameters != electron_mobility_parameters:
                self.set_electron_mobility_parameters(electron_mobility_parameters)

    def set_electron_mobility_parameters(self, electron_mobility_parameters):
        assert isinstance(electron_mobility_parameters, dict), 'Mobility parameters must be dictionary'
        assert 'mu_L0' in electron_mobility_parameters.keys(), 'Mobility parameters must contain mu_L0'
        assert 'v_s' in electron_mobility_parameters.keys(), 'Mobility parameters must contain v_s'
        assert 'A' in electron_mobility_parameters.keys(), 'Mobility parameters must contain A'
        assert 'B' in electron_mobility_parameters.keys(), 'Mobility parameters must contain B'
        assert 'alpha' in electron_mobility_parameters.keys(), 'Mobility parameters must contain alpha'
        assert 'beta' in electron_mobility_parameters.keys(), 'Mobility parameters must contain beta'
        assert isinstance(electron_mobility_parameters['mu_L0'], numbers.Number), 'mu_L0 must be number'
        assert isinstance(electron_mobility_parameters['v_s'], numbers.Number), 'v_s must be number'
        assert isinstance(electron_mobility_parameters['A'], numbers.Number), 'A must be number'
        assert isinstance(electron_mobility_parameters['B'], numbers.Number), 'B must be number'
        assert isinstance(electron_mobility_parameters['alpha'], numbers.Number), 'alpha must be number'
        assert isinstance(electron_mobility_parameters['beta'], numbers.Number), 'beta must be number'
        try:
            mobility_parameter = self.parameters['Electron mobility parameters'].children
            for component in mobility_parameter:
                if component.name == 'mu_L0':
                    component.float_value = float(electron_mobility_parameters['mu_L0'])
                elif component.name == 'v_s':
                    component.float_value = float(electron_mobility_parameters['v_s'])
                elif component.name == 'A':
                    component.float_value = float(electron_mobility_parameters['A'])
                elif component.name == 'B':
                    component.float_value = float(electron_mobility_parameters['B'])
                elif component.name == 'alpha':
                    component.float_value = float(electron_mobility_parameters['alpha'])
                elif component.name == 'beta':
                    component.float_value = float(electron_mobility_parameters['beta'])
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_dict_parameter(name='Electron mobility parameters',
                                                                            description='Electron mobility parameters')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='mu_L0',
                                                                   value=float(electron_mobility_parameters['mu_L0']),
                                                                   description='mu_L0',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='v_s',
                                                                   value=float(electron_mobility_parameters['v_s']),
                                                                   description='v_s',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='A',
                                                                   value=float(electron_mobility_parameters['A']),
                                                                   description='A',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='B',
                                                                   value=float(electron_mobility_parameters['B']),
                                                                   description='B',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='alpha',
                                                                   value=float(electron_mobility_parameters['alpha']),
                                                                   description='alpha parameter',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='beta',
                                                                   value=float(electron_mobility_parameters['beta']),
                                                                   description='beta parameter',
                                                                   parent=parameter)
            self.load_create_sample()

    def _read_in_hole_mobility_parameters(self, hole_mobility_parameters):
        try:
            mobility_parameter = self.parameters['Hole mobility parameters'].children
            for component in mobility_parameter:
                if component.name == 'mu_L0':
                    self.hole_mobility_parameters['mu_L0'] = component.float_value
                elif component.name == 'v_s':
                    self.hole_mobility_parameters['v_s'] = component.float_value
                elif component.name == 'A':
                    self.hole_mobility_parameters['A'] = component.float_value
                elif component.name == 'B':
                    self.hole_mobility_parameters['B'] = component.float_value
                elif component.name == 'alpha':
                    self.hole_mobility_parameters['alpha'] = component.float_value
                elif component.name == 'beta':
                    self.hole_mobility_parameters['beta'] = component.float_value
        except KeyError:
            pass
        if hole_mobility_parameters is not None:
            if self.hole_mobility_parameters != hole_mobility_parameters:
                self.set_hole_mobility_parameters(hole_mobility_parameters)

    def set_hole_mobility_parameters(self, hole_mobility_parameters):
        assert isinstance(hole_mobility_parameters, dict), 'Mobility parameters must be dictionary'
        assert 'mu_L0' in hole_mobility_parameters.keys(), 'Mobility parameters must contain mu_L0'
        assert 'v_s' in hole_mobility_parameters.keys(), 'Mobility parameters must contain v_s'
        assert 'A' in hole_mobility_parameters.keys(), 'Mobility parameters must contain A'
        assert 'B' in hole_mobility_parameters.keys(), 'Mobility parameters must contain B'
        assert 'alpha' in hole_mobility_parameters.keys(), 'Mobility parameters must contain alpha'
        assert 'beta' in hole_mobility_parameters.keys(), 'Mobility parameters must contain beta'
        assert isinstance(hole_mobility_parameters['mu_L0'], numbers.Number), 'mu_L0 must be number'
        assert isinstance(hole_mobility_parameters['v_s'], numbers.Number), 'v_s must be number'
        assert isinstance(hole_mobility_parameters['A'], numbers.Number), 'A must be number'
        assert isinstance(hole_mobility_parameters['B'], numbers.Number), 'B must be number'
        assert isinstance(hole_mobility_parameters['alpha'], numbers.Number), 'alpha must be number'
        assert isinstance(hole_mobility_parameters['beta'], numbers.Number), 'beta must be number'
        try:
            mobility_parameter = self.parameters['Hole mobility parameters'].children
            for component in mobility_parameter:
                if component.name == 'mu_L0':
                    component.float_value = float(hole_mobility_parameters['mu_L0'])
                elif component.name == 'v_s':
                    component.float_value = float(hole_mobility_parameters['v_s'])
                elif component.name == 'A':
                    component.float_value = float(hole_mobility_parameters['A'])
                elif component.name == 'B':
                    component.float_value = float(hole_mobility_parameters['B'])
                elif component.name == 'alpha':
                    component.float_value = float(hole_mobility_parameters['alpha'])
                elif component.name == 'beta':
                    component.float_value = float(hole_mobility_parameters['beta'])
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_dict_parameter(name='Hole mobility parameters',
                                                                            description='Hole mobility parameters')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='mu_L0',
                                                                   value=float(hole_mobility_parameters['mu_L0']),
                                                                   description='mu_L0',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='v_s',
                                                                   value=float(hole_mobility_parameters['v_s']),
                                                                   description='v_s',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='A',
                                                                   value=float(hole_mobility_parameters['A']),
                                                                   description='A',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='B',
                                                                   value=float(hole_mobility_parameters['B']),
                                                                   description='B',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='alpha',
                                                                   value=float(hole_mobility_parameters['alpha']),
                                                                   description='alpha parameter',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='beta',
                                                                   value=float(hole_mobility_parameters['beta']),
                                                                   description='beta parameter',
                                                                   parent=parameter)
            self.load_create_sample()

    def _read_in_thermo_emission_coefficient(self, thermo_emission_coefficient):
        try:
            thermo_emission_parameters = self.parameters['Thermo-emission coefficient'].children
            for particle in thermo_emission_parameters:
                if particle.name == 'electron':
                    self.thermo_emission_coefficient['electron'] = particle.float_value
                elif particle.name == 'hole':
                    self.thermo_emission_coefficient['hole'] = particle.float_value
        except KeyError:
            pass
        if thermo_emission_coefficient is not None:
            if self.thermo_emission_coefficient != thermo_emission_coefficient:
                self.set_thermo_emission_coefficient(thermo_emission_coefficient)

    def set_thermo_emission_coefficient(self, thermo_emission_coefficient):
        assert isinstance(thermo_emission_coefficient, dict), 'Thermo-emission coefficient must be dictionary'
        assert 'electron' in thermo_emission_coefficient.keys(), 'Need values for electron and hole'
        assert 'hole' in thermo_emission_coefficient.keys(), 'Need values for electron and hole'
        assert isinstance(thermo_emission_coefficient['electron'], numbers.Number), 'Coefficient must be number'
        assert isinstance(thermo_emission_coefficient['hole'], numbers.Number), 'Coefficient must be number'
        try:
            thermo_emission_parameters = self.parameters['Thermo-emission coefficient'].children
            for particle in thermo_emission_parameters:
                if particle.name == 'electron':
                    particle.float_value = float(thermo_emission_coefficient['electron'])
                elif particle.name == 'hole':
                    particle.float_value = float(thermo_emission_coefficient['hole'])
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_dict_parameter(name='Thermo-emission coefficient',
                                                                            description='Thermo-emission coefficients')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='electron',
                                                                   value=float(thermo_emission_coefficient['electron']),
                                                                   description='Electron thermo-emission coefficient',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='hole',
                                                                   value=float(thermo_emission_coefficient['hole']),
                                                                   description='Hole thermo-emission coefficient',
                                                                   parent=parameter)
            self.load_create_sample()

    def _read_in_dopants(self, dopants):
        try:
            dopants_list = self.parameters['Dopants'].children
            for dopant in dopants_list:
                dopant_module_name, dopant_class_name = dopant.string_value.split('::')
                dopant_name = dopant.name
                dopant_id = int(dopant.float_value)
                dopant_module = __import__(dopant_module_name, fromlist=[str(dopant_class_name)])
                dopant_class = getattr(dopant_module, dopant_class_name)
                dopant_sample = dopant_class(client=self.client.session_manager, name=dopant_name)
                if dopant_sample.sample.id == dopant_id:
                    self.dopants.append(dopant_sample)
                else:
                    print('Dopant IDs do not match')
        except KeyError:
            pass
        if dopants is not None:
            if self.dopants != dopants:
                self.set_dopants(dopants)

    def set_dopants(self, dopants):
        assert isinstance(dopants, (list, tuple)), 'Expected a list of dopants'
        for dopant in dopants:
            assert isinstance(dopant, Dopant), 'Expected a list of dopants'
            assert isinstance(dopant, Sample), 'Expected a list of dopants'
        try:
            parameter = self.parameters['Dopants']
            for dopant in dopants:
                matched = False
                for dopant_parameter in parameter.children:
                    if dopant.sample.id == int(dopant_parameter.float_value):
                        matched = True
                        dopant_parameter.name = dopant.name
                        dopant_parameter.string_value = dopant.__class__.__module__ + '::' + dopant.__class__.__name__
                        dopant_parameter.description = dopant.description
                        break
                if not matched:
                    self.client.parameter_manager.create_generic_parameter(
                        name=dopant.name,
                        float_value=float(dopant.sample.id),
                        string_value=dopant.__class__.__module__ + '::' + dopant.__class__.__name__,
                        description=dopant.description,
                        parent=parameter)
            dopant_sample_ids = [dopant.sample.id for dopant in dopants]
            for dopant_parameter in parameter.children:
                if int(dopant_parameter.float_value) not in dopant_sample_ids:
                    self.client.parameter_manager.delete_parameter(dopant_parameter)
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_dict_parameter(name='Dopants',
                                                                            description='Dopants dictionary')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            for dopant in dopants:
                self.client.parameter_manager.create_generic_parameter(
                    name=dopant.name,
                    float_value=float(dopant.sample.id),
                    string_value=dopant.__class__.__module__ + '::' + dopant.__class__.__name__,
                    description=dopant.description,
                    parent=parameter)
        self.parameters = {}
        self.load_create_sample()
        self.dopants = dopants

    def _read_in_interfaces(self, interfaces):
        try:
            interfaces_list = self.parameters['Interfaces'].children
            for interface in interfaces_list:
                interface_module_name, interface_class_name = interface.string_value.split('::')
                interface_name = interface.name
                interface_id = int(interface.float_value)
                interface_module = __import__(interface_module_name, fromlist=[interface_class_name])
                interface_class = getattr(interface_module, interface_class_name)
                interface_sample = interface_class(client=self.client.session_manager, name=interface_name)
                if interface_sample.sample.id == interface_id:
                    self.interfaces.append(interface_sample)
                else:
                    print('Interface IDs do not match')
        except KeyError:
            pass
        if interfaces is not None:
            if self.interfaces != interfaces:
                self.set_interfaces(interfaces)

    def set_interfaces(self, interfaces):
        assert isinstance(interfaces, (list, tuple)), 'Expected a list of interfaces'
        for interface in interfaces:
            assert isinstance(interface, Interface), 'Expected a list of interfaces'
            assert isinstance(interface, Sample), 'Expected a list of interfaces'
        try:
            parameter = self.parameters['Interfaces']
            for interface in interfaces:
                matched = False
                for interface_parameter in parameter.children:
                    if interface.sample.id == int(interface_parameter.float_value):
                        matched = True
                        interface_parameter.name = interface.name
                        interface_parameter.string_value = interface.__class__.__module__
                        interface_parameter.string_value += '::' + interface.__class__.__name__
                        interface_parameter.description = interface.description
                        break
                if not matched:
                    self.client.parameter_manager.create_generic_parameter(
                        name=interface.name,
                        float_value=float(interface.sample.id),
                        string_value=interface.__class__.__module__ + '::' + interface.__class__.__name__,
                        description=interface.description,
                        parent=parameter)
            interface_sample_ids = [interface.sample.id for interface in interfaces]
            for interface_parameter in parameter.children:
                if int(interface_parameter.float_value) not in interface_sample_ids:
                    self.client.parameter_manager.delete_parameter(interface_parameter)
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_dict_parameter(name='Interfaces',
                                                                            description='Interfaces dictionary')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            for interface in interfaces:
                self.client.parameter_manager.create_generic_parameter(
                    name=interface.name,
                    float_value=float(interface.sample.id),
                    string_value=interface.__class__.__module__ + '::' + interface.__class__.__name__,
                    description=interface.description,
                    parent=parameter)
        self.parameters = {}
        self.load_create_sample()
        self.interfaces = interfaces
