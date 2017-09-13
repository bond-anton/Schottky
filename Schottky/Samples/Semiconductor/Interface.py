from __future__ import division, print_function
import numpy as np

from Schottky.Samples import Sample


class Interface(Sample):

    def __init__(self, client, name,
                 depth=None, smooth_dirac_epsilon=None,
                 traps=None, description=None):
        super(Interface, self).__init__(client=client, name=name, description=description)
        self.load_create_sample()
        self.__depth = None
        self.__smooth_dirac_epsilon = None
        self.traps = []
        self._read_in_depth(depth)
        self._read_in_smooth_dirac_epsilon(smooth_dirac_epsilon)
        self._read_in_traps(traps)

    @property
    def depth(self):
        return self.__depth

    @depth.setter
    def depth(self, depth):
        try:
            self.parameters['depth'].float_value = np.float64(depth)
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_numeric_parameter(name='depth',
                                                                               value=np.float64(depth),
                                                                               description='Interface depth')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.reload_parameters()
        self.__depth = np.float64(depth)

    def _read_in_depth(self, depth):
        try:
            self.depth = self.parameters['depth'].float_value
        except KeyError:
            pass
        if self.depth != depth and depth is not None:
            self.depth = depth

    @property
    def smooth_dirac_epsilon(self):
        return self.__smooth_dirac_epsilon

    @smooth_dirac_epsilon.setter
    def smooth_dirac_epsilon(self, smooth_dirac_epsilon):
        try:
            self.parameters['smooth dirac epsilon'].float_value = np.float64(smooth_dirac_epsilon)
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_numeric_parameter(name='smooth dirac epsilon',
                                                                               value=np.float64(smooth_dirac_epsilon),
                                                                               description='Smooth dirac epsilon')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.reload_parameters()
        self.__smooth_dirac_epsilon = np.float64(smooth_dirac_epsilon)

    def _read_in_smooth_dirac_epsilon(self, smooth_dirac_epsilon):
        try:
            self.smooth_dirac_epsilon = self.parameters['smooth dirac epsilon'].float_value
        except KeyError:
            pass
        if self.smooth_dirac_epsilon != smooth_dirac_epsilon and smooth_dirac_epsilon is not None:
            self.smooth_dirac_epsilon = smooth_dirac_epsilon

    def _read_in_traps(self, traps):
        pass
