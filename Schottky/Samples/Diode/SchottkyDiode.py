from __future__ import division, print_function
import numpy as np

from Schottky.Samples import Sample


class SchottkyDiode(Sample):

    def __init__(self, client, name,
                 area=None,
                 thickness=None,
                 serial_resistance=None,
                 metal=None,
                 semiconductor=None,
                 description=None):
        super(SchottkyDiode, self).__init__(client=client, name=name, description=description)
        self.load_create_sample()
        self.__area = None
        self.__thickness = None
        self.__serial_resistance = None
        self.__metal = None
        self.__semiconductor = None
        self._read_in_area(area)
        self._read_in_thickness(thickness)
        self._read_in_serial_resistance(serial_resistance)
        self._read_in_metal(metal)
        self._read_in_semiconductor(semiconductor)

    @property
    def area(self):
        return self.__area

    @area.setter
    def area(self, area):
        try:
            self.parameters['area'].float_value = np.float64(area)
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_numeric_parameter(name='area',
                                                                               value=np.float64(area),
                                                                               description='Diode area')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.reload_parameters()
        self.__area = np.float64(area)

    def _read_in_area(self, area):
        try:
            self.area = self.parameters['area'].float_value
        except KeyError:
            pass
        if self.area != area and area is not None:
            self.area = area

    @property
    def thickness(self):
        return self.__thickness

    @thickness.setter
    def thickness(self, thickness):
        try:
            self.parameters['thickness'].float_value = np.float64(thickness)
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_numeric_parameter(name='thickness',
                                                                               value=np.float64(thickness),
                                                                               description='Diode thickness')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.reload_parameters()
        self.__thickness = np.float64(thickness)

    def _read_in_thickness(self, thickness):
        try:
            self.thickness = self.parameters['thickness'].float_value
        except KeyError:
            pass
        if self.thickness != thickness and thickness is not None:
            self.thickness = thickness

    @property
    def serial_resistance(self):
        return self.__serial_resistance

    @serial_resistance.setter
    def serial_resistance(self, serial_resistance):
        try:
            self.parameters['serial resistance'].float_value = np.float64(serial_resistance)
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_numeric_parameter(name='serial resistance',
                                                                               value=np.float64(serial_resistance),
                                                                               description='Diode serial resistance')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.reload_parameters()
        self.__serial_resistance = np.float64(serial_resistance)

    def _read_in_serial_resistance(self, serial_resistance):
        try:
            self.serial_resistance = self.parameters['serial resistance'].float_value
        except KeyError:
            pass
        if self.serial_resistance != serial_resistance and serial_resistance is not None:
            self.serial_resistance = serial_resistance

    @property
    def metal(self):
        return self.__metal

    @metal.setter
    def metal(self, metal):
        string_value = metal.__class__.__module__ + '::'
        string_value += metal.__class__.__name__ + '::'
        string_value += metal.name
        try:
            metal_parameter = self.parameters['Metal']
            metal_parameter.float_value = metal.sample.id
            metal_parameter.string_value = string_value
            metal_parameter.description = metal.description
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_generic_parameter(
                name='Metal',
                float_value=np.float64(metal.sample.id),
                string_value=string_value,
                description=metal.description)
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
        self.parameters = {}
        self.reload_parameters()
        self.__metal = metal

    def _read_in_metal(self, metal):
        try:
            metal_parameter = self.parameters['Metal']
            metal_module_name, metal_class_name, metal_name = metal_parameter.string_value.split('::')
            metal_id = int(metal_parameter.float_value)
            metal_module = __import__(metal_module_name, fromlist=[str(metal_class_name)])
            metal_class = getattr(metal_module, metal_class_name)
            metal_sample = metal_class(client=self.client.session_manager, name=metal_name)
            if metal_sample.sample.id == metal_id:
                self.metal = metal_sample
            else:
                print('Metal IDs do not match')
        except KeyError:
            pass
        if metal is not None:
            if self.metal != metal:
                self.metal = metal

    @property
    def semiconductor(self):
        return self.__semiconductor

    @semiconductor.setter
    def semiconductor(self, semiconductor):
        string_value = semiconductor.__class__.__module__ + '::'
        string_value += semiconductor.__class__.__name__ + '::'
        string_value += semiconductor.name
        try:
            semiconductor_parameter = self.parameters['Semiconductor']
            semiconductor_parameter.float_value = semiconductor.sample.id
            semiconductor_parameter.string_value = string_value
            semiconductor_parameter.description = semiconductor.description
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_generic_parameter(
                name='Semiconductor',
                float_value=np.float64(semiconductor.sample.id),
                string_value=string_value,
                description=semiconductor.description)
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
        self.parameters = {}
        self.reload_parameters()
        self.__semiconductor = semiconductor

    def _read_in_semiconductor(self, semiconductor):
        try:
            semiconductor_parameter = self.parameters['Semiconductor']
            module_name, class_name, semiconductor_name = semiconductor_parameter.string_value.split('::')
            semiconductor_id = int(semiconductor_parameter.float_value)
            semiconductor_module = __import__(module_name, fromlist=[str(class_name)])
            semiconductor_class = getattr(semiconductor_module, class_name)
            semiconductor_sample = semiconductor_class(client=self.client.session_manager, name=semiconductor_name)
            if semiconductor_sample.sample.id == semiconductor_id:
                self.semiconductor = semiconductor_sample
            else:
                print('Semiconductor IDs do not match')
        except KeyError:
            pass
        if semiconductor is not None:
            if self.semiconductor != semiconductor:
                self.semiconductor = semiconductor
