from ScientificProjects.Client import Client
from ScientificProjects.Entities.Sample import Sample as ProjectSample


class Sample(ProjectSample):

    def __init__(self, client, name, description=None):
        assert isinstance(client, Client), 'Valid ScientificProjects Client instance is required'
        self.client = client.user_manager
        super(Sample, self).__init__(name=name,
                                     description=description)
