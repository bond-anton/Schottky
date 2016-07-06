from __future__ import division

from Schottky import constants
from Schottky.Samples import Sample


class Metal(Sample):

    def __init__(self, client, name, work_function=None, description=None):
        super(Metal, self).__init__(client=client, name=name, description=description)
        self.load_create_sample()
        if self.parameters:
            self.work_function = self.parameters['Work function'].float_value
        else:
            self.work_function = None
        if self.work_function != work_function:
            self.set_work_function(work_function)

    def set_work_function(self, work_function):
        if self.parameters:
            self.parameters['Work function'].float_value = work_function
            self.save_sample_changes()
        else:
            parameter = self.client.parameter_manager.create_numeric_parameter(name='Work function',
                                                                               value=work_function,
                                                                               unit_name='eV',
                                                                               description='Metal work function')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.load_create_sample()
        self.work_function = work_function
