from __future__ import division, print_function
import numbers
import numpy as np

from Schottky.Samples import Sample
from Schottky.Samples.Semiconductor.Dopant import Dopant
from Schottky.Samples.Semiconductor.Interface import Interface


class Semiconductor(Sample):
    def __init__(self, client, name,
                 epsilon=None,
                 affinity=None,
                 effective_mass=None,
                 effective_bands_density_of_states=None,
                 band_gap_parameters=None,
                 electron_mobility_parameters=None,
                 hole_mobility_parameters=None,
                 thermo_emission_coefficient=None,
                 dopants=None, interfaces=None,
                 description=None):
        super(Semiconductor, self).__init__(client=client, name=name, description=description)
        self.load_create_sample()
        self.__epsilon = None
        self.__affinity = None
        self.__effective_mass = {}
        self.__effective_bands_density_of_states = {}
        self.__band_gap_parameters = {}
        self.__electron_mobility_parameters = {}
        self.__hole_mobility_parameters = {}
        self.__thermo_emission_coefficient = {}
        self.__dopants = []
        self.__interfaces = []
        self._read_in_epsilon(epsilon)
        self._read_in_affinity(affinity)
        self._read_in_effective_mass(effective_mass)
        self._read_in_effective_bands_density_of_states(effective_bands_density_of_states)
        self._read_in_band_gap_parameters(band_gap_parameters)
        self._read_in_electron_mobility_parameters(electron_mobility_parameters)
        self._read_in_hole_mobility_parameters(hole_mobility_parameters)
        self._read_in_thermo_emission_coefficient(thermo_emission_coefficient)
        self._read_in_dopants(dopants)
        self._read_in_interfaces(interfaces)

    @property
    def epsilon(self):
        return self.__epsilon

    @epsilon.setter
    def epsilon(self, epsilon):
        try:
            self.parameters['epsilon'].float_value = np.float64(epsilon)
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_numeric_parameter(name='epsilon',
                                                                               value=np.float64(epsilon),
                                                                               description='Permittivity of space')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.reload_parameters()
        self.__epsilon = np.float64(epsilon)

    def _read_in_epsilon(self, epsilon):
        try:
            self.__epsilon = np.float64(self.parameters['epsilon'].float_value)
        except KeyError:
            pass
        if self.epsilon != epsilon and epsilon is not None:
            self.epsilon = epsilon

    @property
    def affinity(self):
        return self.__affinity

    @affinity.setter
    def affinity(self, affinity):
        try:
            self.parameters['affinity'].float_value = np.float64(affinity)
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_numeric_parameter(name='affinity',
                                                                               value=np.float64(affinity),
                                                                               description='Electron affinity')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.reload_parameters()
        self.__affinity = np.float64(affinity)

    def _read_in_affinity(self, affinity):
        try:
            self.__affinity = np.float64(self.parameters['affinity'].float_value)
        except KeyError:
            pass
        if self.affinity != affinity and affinity is not None:
            self.affinity = affinity

    @property
    def effective_mass(self):
        return self.__effective_mass

    @effective_mass.setter
    def effective_mass(self, effective_mass):
        try:
            effective_mass_parameters = self.parameters['Effective mass'].children
            for particle in effective_mass_parameters:
                if particle.name == 'electron':
                    particle.float_value = np.float64(effective_mass['electron'])
                elif particle.name == 'hole':
                    particle.float_value = np.float64(effective_mass['hole'])
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_dict_parameter(name='Effective mass',
                                                                            description='Effective masses of carriers')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='electron',
                                                                   value=np.float64(effective_mass['electron']),
                                                                   description='Electron effective mass',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='hole',
                                                                   value=np.float64(effective_mass['hole']),
                                                                   description='Hole effective mass',
                                                                   parent=parameter)
            self.reload_parameters()
        self.__effective_mass = {'electron': np.float64(effective_mass['electron']),
                                 'hole': np.float64(effective_mass['hole'])}

    def _read_in_effective_mass(self, effective_mass):
        try:
            effective_mass_parameters = self.parameters['Effective mass'].children
            for particle in effective_mass_parameters:
                if particle.name == 'electron':
                    self.__effective_mass['electron'] = particle.float_value
                elif particle.name == 'hole':
                    self.__effective_mass['hole'] = particle.float_value
        except KeyError:
            pass
        if effective_mass is not None:
            if self.effective_mass != effective_mass:
                self.effective_mass = effective_mass

    @property
    def effective_bands_density_of_states(self):
        return self.__effective_bands_density_of_states

    @effective_bands_density_of_states.setter
    def effective_bands_density_of_states(self, effective_bands_density_of_states):
        try:
            bands_parameters = self.parameters['Bands density of states'].children
            for band in bands_parameters:
                if band.name == 'Nc':
                    band.float_value = np.float64(effective_bands_density_of_states['Nc'])
                elif band.name == 'Nv':
                    band.float_value = np.float64(effective_bands_density_of_states['Nv'])
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_dict_parameter(name='Bands density of states',
                                                                            description='Bands density of states')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='Nc',
                                                                   value=np.float64(
                                                                       effective_bands_density_of_states['Nc']),
                                                                   description='Nc density of states',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='Nv',
                                                                   value=np.float64(
                                                                       effective_bands_density_of_states['Nv']),
                                                                   description='Nv density of states',
                                                                   parent=parameter)
            self.reload_parameters()
            self.__effective_bands_density_of_states = {'Nc': np.float64(effective_bands_density_of_states['Nc']),
                                                        'Nv': np.float64(effective_bands_density_of_states['Nv'])}

    def _read_in_effective_bands_density_of_states(self, effective_bands_density_of_states):
        try:
            bands_parameters = self.parameters['Bands density of states'].children
            for band in bands_parameters:
                if band.name == 'Nc':
                    self.__effective_bands_density_of_states['Nc'] = np.float64(band.float_value)
                elif band.name == 'Nv':
                    self.__effective_bands_density_of_states['Nv'] = np.float64(band.float_value)
        except KeyError:
            pass
        if effective_bands_density_of_states is not None:
            if self.effective_bands_density_of_states != effective_bands_density_of_states:
                self.effective_bands_density_of_states = effective_bands_density_of_states

    @property
    def band_gap_parameters(self):
        return self.__band_gap_parameters

    @band_gap_parameters.setter
    def band_gap_parameters(self, band_gap_parameters):
        try:
            band_gap_parameter = self.parameters['Band gap parameters'].children
            for component in band_gap_parameter:
                if component.name == 'Eg_0':
                    component.float_value = np.float64(band_gap_parameters['Eg_0'])
                elif component.name == 'alpha':
                    component.float_value = np.float64(band_gap_parameters['alpha'])
                elif component.name == 'beta':
                    component.float_value = np.float64(band_gap_parameters['beta'])
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_dict_parameter(name='Band gap parameters',
                                                                            description='Band gap parameters')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='Eg_0',
                                                                   value=np.float64(band_gap_parameters['Eg_0']),
                                                                   description='Eg_0',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='alpha',
                                                                   value=np.float64(band_gap_parameters['alpha']),
                                                                   description='alpha parameter',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='beta',
                                                                   value=np.float64(band_gap_parameters['beta']),
                                                                   description='beta parameter',
                                                                   parent=parameter)
            self.reload_parameters()
            self.__band_gap_parameters = {'Eg_0': np.float64(band_gap_parameters['Eg_0']),
                                          'alpha': np.float64(band_gap_parameters['alpha']),
                                          'beta': np.float64(band_gap_parameters['beta'])}

    def _read_in_band_gap_parameters(self, band_gap_parameters):
        try:
            band_gap_parameter = self.parameters['Band gap parameters'].children
            for component in band_gap_parameter:
                if component.name == 'Eg_0':
                    self.__band_gap_parameters['Eg_0'] = np.float64(component.float_value)
                elif component.name == 'alpha':
                    self.__band_gap_parameters['alpha'] = np.float64(component.float_value)
                elif component.name == 'beta':
                    self.__band_gap_parameters['beta'] = np.float64(component.float_value)
        except KeyError:
            pass
        if band_gap_parameters is not None:
            if self.band_gap_parameters != band_gap_parameters:
                self.band_gap_parameters = band_gap_parameters

    @property
    def electron_mobility_parameters(self):
        return self.__electron_mobility_parameters

    @electron_mobility_parameters.setter
    def electron_mobility_parameters(self, electron_mobility_parameters):
        try:
            mobility_parameter = self.parameters['Electron mobility parameters'].children
            for component in mobility_parameter:
                if component.name == 'mu_L0':
                    component.float_value = np.float64(electron_mobility_parameters['mu_L0'])
                elif component.name == 'v_s':
                    component.float_value = np.float64(electron_mobility_parameters['v_s'])
                elif component.name == 'A':
                    component.float_value = np.float64(electron_mobility_parameters['A'])
                elif component.name == 'B':
                    component.float_value = np.float64(electron_mobility_parameters['B'])
                elif component.name == 'alpha':
                    component.float_value = np.float64(electron_mobility_parameters['alpha'])
                elif component.name == 'beta':
                    component.float_value = np.float64(electron_mobility_parameters['beta'])
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_dict_parameter(name='Electron mobility parameters',
                                                                            description='Electron mobility parameters')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='mu_L0',
                                                                   value=np.float64(
                                                                       electron_mobility_parameters['mu_L0']),
                                                                   description='mu_L0',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='v_s',
                                                                   value=np.float64(
                                                                       electron_mobility_parameters['v_s']),
                                                                   description='v_s',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='A',
                                                                   value=np.float64(electron_mobility_parameters['A']),
                                                                   description='A',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='B',
                                                                   value=np.float64(electron_mobility_parameters['B']),
                                                                   description='B',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='alpha',
                                                                   value=np.float64(
                                                                       electron_mobility_parameters['alpha']),
                                                                   description='alpha parameter',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='beta',
                                                                   value=np.float64(
                                                                       electron_mobility_parameters['beta']),
                                                                   description='beta parameter',
                                                                   parent=parameter)
            self.reload_parameters()
            self.__electron_mobility_parameters = {'mu_L0': np.float64(electron_mobility_parameters['mu_L0']),
                                                   'v_s': np.float64(electron_mobility_parameters['v_s']),
                                                   'A': np.float64(electron_mobility_parameters['A']),
                                                   'B': np.float64(electron_mobility_parameters['B']),
                                                   'beta': np.float64(electron_mobility_parameters['beta'])}

    def _read_in_electron_mobility_parameters(self, electron_mobility_parameters):
        try:
            mobility_parameter = self.parameters['Electron mobility parameters'].children
            for component in mobility_parameter:
                if component.name == 'mu_L0':
                    self.__electron_mobility_parameters['mu_L0'] = np.float64(component.float_value)
                elif component.name == 'v_s':
                    self.__electron_mobility_parameters['v_s'] = np.float64(component.float_value)
                elif component.name == 'A':
                    self.__electron_mobility_parameters['A'] = np.float64(component.float_value)
                elif component.name == 'B':
                    self.__electron_mobility_parameters['B'] = np.float64(component.float_value)
                elif component.name == 'alpha':
                    self.__electron_mobility_parameters['alpha'] = np.float64(component.float_value)
                elif component.name == 'beta':
                    self.__electron_mobility_parameters['beta'] = np.float64(component.float_value)
        except KeyError:
            pass
        if electron_mobility_parameters is not None:
            if self.electron_mobility_parameters != electron_mobility_parameters:
                self.electron_mobility_parameters = electron_mobility_parameters

    @property
    def hole_mobility_parameters(self):
        return self.__hole_mobility_parameters

    @hole_mobility_parameters.setter
    def hole_mobility_parameters(self, hole_mobility_parameters):
        try:
            mobility_parameter = self.parameters['Hole mobility parameters'].children
            for component in mobility_parameter:
                if component.name == 'mu_L0':
                    component.float_value = np.float64(hole_mobility_parameters['mu_L0'])
                elif component.name == 'v_s':
                    component.float_value = np.float64(hole_mobility_parameters['v_s'])
                elif component.name == 'A':
                    component.float_value = np.float64(hole_mobility_parameters['A'])
                elif component.name == 'B':
                    component.float_value = np.float64(hole_mobility_parameters['B'])
                elif component.name == 'alpha':
                    component.float_value = np.float64(hole_mobility_parameters['alpha'])
                elif component.name == 'beta':
                    component.float_value = np.float64(hole_mobility_parameters['beta'])
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_dict_parameter(name='Hole mobility parameters',
                                                                            description='Hole mobility parameters')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='mu_L0',
                                                                   value=np.float64(hole_mobility_parameters['mu_L0']),
                                                                   description='mu_L0',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='v_s',
                                                                   value=np.float64(hole_mobility_parameters['v_s']),
                                                                   description='v_s',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='A',
                                                                   value=np.float64(hole_mobility_parameters['A']),
                                                                   description='A',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='B',
                                                                   value=np.float64(hole_mobility_parameters['B']),
                                                                   description='B',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='alpha',
                                                                   value=np.float64(hole_mobility_parameters['alpha']),
                                                                   description='alpha parameter',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='beta',
                                                                   value=np.float64(hole_mobility_parameters['beta']),
                                                                   description='beta parameter',
                                                                   parent=parameter)
            self.reload_parameters()
            self.__hole_mobility_parameters = {'mu_L0': np.float64(hole_mobility_parameters['mu_L0']),
                                               'v_s': np.float64(hole_mobility_parameters['v_s']),
                                               'A': np.float64(hole_mobility_parameters['A']),
                                               'B': np.float64(hole_mobility_parameters['B']),
                                               'alpha': np.float64(hole_mobility_parameters['alpha']),
                                               'beta': np.float64(hole_mobility_parameters['beta'])}

    def _read_in_hole_mobility_parameters(self, hole_mobility_parameters):
        try:
            mobility_parameter = self.parameters['Hole mobility parameters'].children
            for component in mobility_parameter:
                if component.name == 'mu_L0':
                    self.__hole_mobility_parameters['mu_L0'] = np.float64(component.float_value)
                elif component.name == 'v_s':
                    self.__hole_mobility_parameters['v_s'] = np.float64(component.float_value)
                elif component.name == 'A':
                    self.__hole_mobility_parameters['A'] = np.float64(component.float_value)
                elif component.name == 'B':
                    self.__hole_mobility_parameters['B'] = np.float64(component.float_value)
                elif component.name == 'alpha':
                    self.__hole_mobility_parameters['alpha'] = np.float64(component.float_value)
                elif component.name == 'beta':
                    self.__hole_mobility_parameters['beta'] = np.float64(component.float_value)
        except KeyError:
            pass
        if hole_mobility_parameters is not None:
            if self.hole_mobility_parameters != hole_mobility_parameters:
                self.hole_mobility_parameters = hole_mobility_parameters

    @property
    def thermo_emission_coefficient(self):
        return self.__thermo_emission_coefficient

    @thermo_emission_coefficient.setter
    def thermo_emission_coefficient(self, thermo_emission_coefficient):
        try:
            thermo_emission_parameters = self.parameters['Thermo-emission coefficient'].children
            for particle in thermo_emission_parameters:
                if particle.name == 'electron':
                    particle.float_value = np.float64(thermo_emission_coefficient['electron'])
                elif particle.name == 'hole':
                    particle.float_value = np.float64(thermo_emission_coefficient['hole'])
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_dict_parameter(name='Thermo-emission coefficient',
                                                                            description='Thermo-emission coefficients')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='electron',
                                                                   value=np.float64(thermo_emission_coefficient['electron']),
                                                                   description='Electron thermo-emission coefficient',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='hole',
                                                                   value=np.float64(thermo_emission_coefficient['hole']),
                                                                   description='Hole thermo-emission coefficient',
                                                                   parent=parameter)
            self.reload_parameters()
            self.__thermo_emission_coefficient = {'electron': np.float64(thermo_emission_coefficient['electron']),
                                                  'hole': np.float64(thermo_emission_coefficient['hole'])}

    def _read_in_thermo_emission_coefficient(self, thermo_emission_coefficient):
        try:
            thermo_emission_parameters = self.parameters['Thermo-emission coefficient'].children
            for particle in thermo_emission_parameters:
                if particle.name == 'electron':
                    self.__thermo_emission_coefficient['electron'] = particle.float_value
                elif particle.name == 'hole':
                    self.__thermo_emission_coefficient['hole'] = particle.float_value
        except KeyError:
            pass
        if thermo_emission_coefficient is not None:
            if self.thermo_emission_coefficient != thermo_emission_coefficient:
                self.thermo_emission_coefficient = thermo_emission_coefficient

    @property
    def dopants(self):
        return self.__dopants

    @dopants.setter
    def dopants(self, dopants):
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
                        float_value=np.float64(dopant.sample.id),
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
                    float_value=np.float64(dopant.sample.id),
                    string_value=dopant.__class__.__module__ + '::' + dopant.__class__.__name__,
                    description=dopant.description,
                    parent=parameter)
        self.parameters = {}
        self.reload_parameters()
        self.__dopants = dopants

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
                    self.__dopants.append(dopant_sample)
                else:
                    print('Dopant IDs do not match')
        except KeyError:
            pass
        if dopants is not None:
            if self.dopants != dopants:
                self.dopants = dopants

    @property
    def interfaces(self):
        return self.__interfaces

    @interfaces.setter
    def interfaces(self, interfaces):
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
                        float_value=np.float64(interface.sample.id),
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
                    float_value=np.float64(interface.sample.id),
                    string_value=interface.__class__.__module__ + '::' + interface.__class__.__name__,
                    description=interface.description,
                    parent=parameter)
        self.parameters = {}
        self.reload_parameters()
        self.__interfaces = interfaces

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
                    self.__interfaces.append(interface_sample)
                else:
                    print('Interface IDs do not match')
        except KeyError:
            pass
        if interfaces is not None:
            if self.interfaces != interfaces:
                self.interfaces = interfaces
