from __future__ import division, print_function
import numpy as np

from ScientificProjects.Client import Client
from ScientificProjects.Entities.Parameter import Parameter


class Sample(object):

    def __init__(self, client, name, description=None, parameters=None):
        assert isinstance(client, Client), 'Valid ScientificProjects Client instance is required'
        assert client.user_manager.project_manager.project_opened(), 'Please open project to work with.'
        self.client = client.user_manager
        self.name = name
        self.description = description
        self.parameters = {}
        if isinstance(parameters, (list, tuple, np.ndarray)):
            for parameter in parameters:
                if isinstance(parameter, Parameter):
                    self.parameters[parameter.name] = parameter
        self.sample = None

    def load_create_sample(self):
        self.sample = self.client.sample_manager.create_sample(name=self.name, description=self.description)
        if self.sample.parameters:
            for parameter in self.sample.parameters:
                try:
                    if self.parameters[parameter.name] != parameter:
                        raise ValueError('To change sample parameter first open it without providing any parameters')
                except KeyError:
                    self.parameters[parameter.name] = parameter
            for parameter in self.parameters.values():
                if parameter not in self.sample.parameters:
                    raise ValueError('To add new parameter to sample first open it without providing any parameters')
        else:
            for parameter in self.parameters.values():
                self.client.sample_manager.add_parameter_to_sample(self.sample, parameter)

    def save_sample_changes(self):
        self.client.session.commit()

    def __str__(self):
        description = str(self.sample)
        for parameter in self.parameters.values():
            description += '\n' + str(parameter)
        return description
