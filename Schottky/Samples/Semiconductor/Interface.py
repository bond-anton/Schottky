from __future__ import division, print_function
import numbers

from Schottky.Samples import Sample
from Schottky.Samples.Trap import Trap


class Interface(Sample):

    def __init__(self, client, name,
                 depth=None, smooth_dirac_epsilon=None,
                 traps=None, description=None):
        super(Interface, self).__init__(client=client, name=name, description=description)
        self.load_create_sample()
        self.depth = None
        self.smooth_dirac_epsilon = None
        self.traps = []
        self._read_in_depth(depth)
        self._read_in_smooth_dirac_epsilon(smooth_dirac_epsilon)
        self._read_in_traps(traps)

    def _read_in_depth(self, depth):
        try:
            self.depth = self.parameters['depth'].float_value
        except KeyError:
            pass
        if self.depth != depth and depth is not None:
            self.set_depth(depth)

    def set_depth(self, depth):
        assert isinstance(depth, numbers.Number), 'Depth must be a number'
        try:
            self.parameters['depth'].float_value = float(depth)
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_numeric_parameter(name='depth',
                                                                               value=float(depth),
                                                                               description='Interface depth')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.load_create_sample()
        self.depth = float(depth)

    def _read_in_smooth_dirac_epsilon(self, smooth_dirac_epsilon):
        try:
            self.depth = self.parameters['smooth dirac epsilon'].float_value
        except KeyError:
            pass
        if self.smooth_dirac_epsilon != smooth_dirac_epsilon and smooth_dirac_epsilon is not None:
            self.set_smooth_dirac_epsilon(smooth_dirac_epsilon)

    def set_smooth_dirac_epsilon(self, smooth_dirac_epsilon):
        assert isinstance(smooth_dirac_epsilon, numbers.Number), 'Smooth dirac epsilon must be a number'
        try:
            self.parameters['smooth dirac epsilon'].float_value = float(smooth_dirac_epsilon)
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_numeric_parameter(name='smooth dirac epsilon',
                                                                               value=float(smooth_dirac_epsilon),
                                                                               description='Smooth dirac epsilon')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.load_create_sample()
        self.smooth_dirac_epsilon = float(smooth_dirac_epsilon)

    def _read_in_traps(self, traps):
        pass
