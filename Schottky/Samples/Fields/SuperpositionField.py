from __future__ import division, print_function

from sqlalchemy.orm.exc import StaleDataError

from Space.Field import Field, SuperposedField

from Schottky.Samples import Sample


class SuperpositionField(Sample, SuperposedField):

    def __init__(self, client, name, fields=None, description=None):
        Sample.__init__(self, client=client, name=name, description=description)
        self.load_create_sample()
        self.fields = []
        self._read_in_fields(fields=fields)
        SuperposedField.__init__(self, name=name, fields=self.fields)

    def _read_in_fields(self, fields):
        try:
            field_components = self.parameters['Fields'].children
            for field in field_components:
                field_module_name, field_class_name = field.string_value.split('::')
                field_name = field.name
                field_id = int(field.float_value)
                field_module = __import__(field_module_name, fromlist=[field_class_name])
                field_class = getattr(field_module, field_class_name)
                field_sample = field_class(client=self.client.session_manager, name=field_name)
                if field_sample.sample.id == field_id:
                    self.fields.append(field_sample)
                else:
                    print('Field IDs do not match')
        except KeyError:
            pass
        if fields is not None:
            if self.fields != fields:
                self.set_fields(fields)

    def set_fields(self, fields):
        assert isinstance(fields, (list, tuple)), 'Expected a list of fields'
        for field in fields:
            assert isinstance(field, Field), 'Expected a list of fields'
            assert isinstance(field, Sample), 'Expected a list of fields'
        try:
            for field in self.parameters['Fields'].children:
                self.client.parameter_manager.delete_parameter(field)
            try:
                self.client.parameter_manager.delete_parameter(self.parameters['Fields'])
            except StaleDataError:
                self.client.session.rollback()
                pass
            self.save_sample_changes()
        except KeyError:
            pass
        parameter = self.client.parameter_manager._create_parameter(name='Fields',
                                                                    description='Fields dictionary',
                                                                    parameter_type='Dictionary')
        self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                           parameter=parameter)
        for field in fields:
            self.client.parameter_manager._create_parameter(
                name=field.name,
                parameter_type='Generic',
                float_value=float(field.sample.id),
                string_value=field.__class__.__module__ + '::' + field.__class__.__name__,
                description=field.description,
                parent=parameter)
        self.parameters = {}
        self.load_create_sample()
        self.fields = fields
