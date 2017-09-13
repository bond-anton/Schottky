from __future__ import division
import numbers
import numpy as np

from Schottky.Samples import Sample


class Metal(Sample):

    def __init__(self, client, name, work_function=None, description=None):
        super(Metal, self).__init__(client=client, name=name, description=description)
        self.load_create_sample()
        self.__work_function = None
        self._read_in_work_function(work_function)

    @property
    def work_function(self):
        return self.__work_function

    @work_function.setter
    def work_function(self, work_function):
        try:
            self.parameters['Work function'].float_value = np.float64(work_function)
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_numeric_parameter(name='Work function',
                                                                               value=work_function,
                                                                               unit_name='eV',
                                                                               description='Metal work function')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.reload_parameters()
        self.__work_function = np.float64(work_function)

    def _read_in_work_function(self, work_function):
        try:
            self.__work_function = np.float64(self.parameters['Work function'].float_value)
        except KeyError:
            pass
        if self.work_function != work_function and work_function is not None:
            self.work_function = work_function
