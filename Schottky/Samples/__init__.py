from __future__ import division, print_function
import numpy as np

from BDProjects.Client import Client
from BDProjects.Entities.Parameter import Parameter


class Sample(object):

    def __init__(self, client, name, description=None, parameters=None):
        assert isinstance(client, Client), 'Valid ScientificProjects Client instance is required'
        assert client.user_manager.project_manager.project_opened(), 'Please open project to work with.'
        self.__client = client.user_manager
        self.__name = str(name)
        self.__description = str(description)
        self.__parameters = {}
        self.parameters = parameters
        self.sample = None

    @property
    def client(self):
        return self.__client

    @property
    def name(self):
        return self.__name

    @name.setter
    def name(self, name):
        self.__name = str(name)

    @property
    def description(self):
        return self.__description

    @description.setter
    def description(self, description):
        self.__description = str(description)

    @property
    def parameters(self):
        return self.__parameters

    @parameters.setter
    def parameters(self, parameters):
        if isinstance(parameters, (list, tuple, np.ndarray)):
            for parameter in parameters:
                if isinstance(parameter, Parameter):
                    self.__parameters[parameter.name] = parameter
        elif isinstance(parameters, dict):
            for parameter_key in parameters.keys():
                if isinstance(parameters[parameter_key], Parameter):
                    self.__parameters[parameters[parameter_key].name] = parameters[parameter_key]
        elif parameters is None:
            pass
        else:
            raise TypeError('Parameters must be provided either as an iterable or as a dictionary')

    def load_create_sample(self):
        samples = self.client.sample_manager.get_samples(name=self.name, exact=True)
        if len(samples) == 1:
            self.sample = samples[0]
        elif len(samples) == 0:
            self.sample = self.client.sample_manager.create_sample(name=self.name, description=self.description)
        else:
            raise ValueError('More than one (%d) sample found for given name, check the database' % len(samples))
        self.reload_parameters()

    def reload_parameters(self):
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
