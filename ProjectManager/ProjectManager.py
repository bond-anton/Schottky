from __future__ import division, print_function
# import time
import datetime
import numpy as np
import mpmath
from scipy.interpolate import interp1d

from sqlalchemy import create_engine, MetaData, Table, Column, Integer, String, Text, DateTime, ForeignKey, BigInteger
from sqlalchemy import desc
from sqlalchemy.sql import func, select, and_, or_

from .numeric_operations import number_to_four_ints


class Project():
    def __init__(self, db_name='/:memory:', backend='sqlite',
                 hostname='', port='', user='', password='',
                 overwrite=False, debug=False):
        if debug: print(db_name)
        if debug: print(backend + '://' + hostname + '/' + db_name)
        # postgresql://scott:tiger@localhost:5432/mydatabase
        credentials = user + ':' + password if password else user
        if credentials:
            credentials += '@'
        if port:
            hostname += ':' + str(port)
        self.engine = create_engine(backend + '://' + credentials + hostname + '/' + db_name)
        self.metadata = MetaData()

        self.Project_Info_Fields = Table('project_info_fields', self.metadata,
                                         Column('field_id', Integer, primary_key=True),
                                         Column('field_name', String(64), nullable=False),
                                         Column('field_description', String(255), nullable=True)
                                         )

        self.Project_Info = Table('project_info', self.metadata,
                                  Column('info_id', Integer, primary_key=True),
                                  Column('info_field_id', Integer, ForeignKey('project_info_fields.field_id')),
                                  Column('info_string', String(255), nullable=False),
                                  Column('info_created_date', DateTime, default=datetime.datetime.now(), nullable=False)
                                  )

        self.Log_Types = Table('log_types', self.metadata,
                               Column('log_type_id', Integer, primary_key=True),
                               Column('log_type_name_long', String(200), nullable=False),
                               Column('log_type_name_short', String(5), nullable=False),
                               Column('log_type_description', Text, nullable=True)
                               )

        self.Log = Table('log', self.metadata,
                         Column('log_id', Integer, primary_key=True),
                         Column('log_type_id', Integer, ForeignKey('log_types.log_type_id')),
                         Column('log_datetime', DateTime, default=datetime.datetime.now(), nullable=False),
                         Column('log_message', Text, nullable=False)
                         )

        self.Tools = Table('tools', self.metadata,
                           Column('tool_id', Integer, primary_key=True),
                           Column('tool_name', String(200), nullable=False),
                           Column('tool_description', Text, nullable=True)
                           )

        self.Measurements_Types = Table('measurements_types', self.metadata,
                                        Column('measurement_type_id', Integer, primary_key=True),
                                        Column('measurement_type_name', String(200), nullable=False),
                                        Column('measurement_type_description', Text, nullable=True)
                                        )

        self.Measurements = Table('measurements', self.metadata,
                                  Column('measurement_id', Integer, primary_key=True),
                                  Column('measurement_type_id', Integer,
                                         ForeignKey('measurements_types.measurement_type_id')),
                                  Column('measurement_name', String(200), nullable=False),
                                  Column('measurement_description', Text, nullable=True),
                                  Column('measurement_datetime', DateTime, nullable=False),
                                  Column('measurement_added', DateTime, default=datetime.datetime.now()),
                                  Column('measurement_altered', DateTime, default=datetime.datetime.now(),
                                         onupdate=datetime.datetime.now())
                                  )

        self.DataPoints = Table('data_points', self.metadata,
                                Column('point_id', Integer, primary_key=True),
                                Column('point_order', Integer, nullable=False),
                                Column('measurement_id', Integer, ForeignKey('measurements.measurement_id')),
                                Column('tool_id', Integer, ForeignKey('tools.tool_id')),
                                Column('point_name_long', String(60), nullable=False),
                                Column('point_name_short', String(20), nullable=False),
                                Column('point_unit_name', String(20), nullable=False),
                                Column('point_sign', Integer, nullable=False),
                                Column('point_mantissa', BigInteger, nullable=False),
                                Column('point_exponent', Integer, nullable=False),
                                Column('point_bytecount', Integer, nullable=False),
                                Column('point_measured', DateTime, nullable=False),
                                Column('point_added', DateTime, default=datetime.datetime.now()),
                                Column('point_altered', DateTime, default=datetime.datetime.now(),
                                       onupdate=datetime.datetime.now())
                                )

        self.create_tables(overwrite)
        self.connection = self.engine.connect()
        self.add_default_log_types()

    def __enter__(self):
        return self

    def __exit__(self):
        self.connection.close()

    def create_tables(self, overwrite=False):
        if overwrite:
            self.metadata.drop_all(self.engine)
        self.metadata.create_all(self.engine)

    def add_project_info_field(self, field_name, field_description=None):
        sel = select([self.Project_Info_Fields]).where(self.Project_Info_Fields.c.field_name.like(field_name))
        result = self.connection.execute(sel)
        for row in result:
            return row['field_id']
        ins = self.Project_Info_Fields.insert().values(field_name=field_name,
                                                       field_description=field_description)
        ins.compile().params
        result = self.connection.execute(ins)
        return result.inserted_primary_key[0]

    def add_project_info(self, field_name, info_string, create_field=True):
        try:
            sel = select([self.Project_Info_Fields]).where(self.Project_Info_Fields.c.field_name.like(field_name))
            result = self.connection.execute(sel)
            row = result.fetchone()
            field_id = row['field_id']
            field_name_rec = row['field_name']
        except:
            if create_field:
                field_id = self.add_project_info_field(field_name)
                field_name_rec = field_name
        ins = self.Project_Info.insert().values(info_string=info_string, info_field_id=field_id)
        ins.compile().params
        result = self.connection.execute(ins)
        self.log_record('Project info field \'%s\' is set to \'%s\'' % (field_name_rec, info_string), 'Information')
        return result.inserted_primary_key

    def add_project_info_if_changed(self, field_name, info_string):
        try:
            info_rec = self.get_project_info(field_name)
        except:
            info_rec = None
        if info_string != info_rec:
            self.add_project_info(field_name, info_string)
            info_rec = info_string
        return info_rec

    def get_project_info(self, field_name):
        sel = select([self.Project_Info_Fields]).where(self.Project_Info_Fields.c.field_name.like(field_name))
        result = self.connection.execute(sel)
        row = result.fetchone()
        field_id = row['field_id']
        # field_name_rec = row['field_name']
        sel = select([self.Project_Info]).where(self.Project_Info.c.info_field_id == field_id).order_by(
            'info_created_date')
        result = self.connection.execute(sel)
        row = result.fetchone()
        return row['info_string']

    def add_log_type(self, log_type_name_long, log_type_name_short, log_type_description=None):
        sel = select([self.Log_Types]).where(self.Log_Types.c.log_type_name_long.like(log_type_name_long))
        result = self.connection.execute(sel)
        for row in result:
            return row['log_type_id']
        ins = self.Log_Types.insert().values(log_type_name_long=log_type_name_long,
                                             log_type_name_short=log_type_name_short,
                                             log_type_description=log_type_description)
        ins.compile().params
        result = self.connection.execute(ins)
        return result.inserted_primary_key[0]

    def add_default_log_types(self):
        self.add_log_type('Information', 'I', 'Information log messages')
        self.add_log_type('Warning', 'W', 'Warning log messages')
        self.add_log_type('Error', 'E', 'Error log messages')
        self.add_log_type('Debug', 'D', 'Debug log messages')

    def log_record(self, message, log_type):
        sel = select([self.Log_Types]).where(self.Log_Types.c.log_type_name_long.like(log_type))
        result = self.connection.execute(sel)
        row = result.fetchone()
        log_type_id = row['log_type_id']
        ins = self.Log.insert().values(log_message=message, log_type_id=log_type_id)
        ins.compile().params
        result = self.connection.execute(ins)
        return result.inserted_primary_key

    def add_tool(self, tool_name, tool_description=None):
        sel = select([self.Tools]).where(self.Tools.c.tool_name == tool_name)
        result = self.connection.execute(sel)
        for row in result:
            return row['tool_id']
        ins = self.Tools.insert().values(tool_name=tool_name, tool_description=tool_description)
        ins.compile().params
        result = self.connection.execute(ins)
        self.log_record('Tool: %s added' % (tool_name), 'Information')
        return result.inserted_primary_key[0]

    def add_measurement_type(self, measurement_type_name, measurement_type_description=None):
        sel = select([self.Measurements_Types]).where(
            self.Measurements_Types.c.measurement_type_name == measurement_type_name)
        result = self.connection.execute(sel)
        for row in result:
            return row['measurement_type_id']
        ins = self.Measurements_Types.insert().values(measurement_type_name=measurement_type_name,
                                                      measurement_type_description=measurement_type_description)
        ins.compile().params
        result = self.connection.execute(ins)
        self.log_record('Measurement type: %s added' % (measurement_type_name), 'Information')
        return result.inserted_primary_key[0]

    def add_measurement(self, measurement_type_id, measurement_name,
                        measurement_description=None, measurement_datetime=None, create_new=True):
        if measurement_datetime is None:
            measurement_datetime = datetime.datetime.now()
        if create_new:
            sel = select([self.Measurements]).where(and_(self.Measurements.c.measurement_name == measurement_name,
                                                         self.Measurements.c.measurement_type_id == measurement_type_id,
                                                         self.Measurements.c.measurement_datetime == measurement_datetime))
        else:
            sel = select([self.Measurements]).where(and_(self.Measurements.c.measurement_name == measurement_name,
                                                         self.Measurements.c.measurement_type_id == measurement_type_id))
        result = self.connection.execute(sel)
        for row in result:
            return row['measurement_id']
        ins = self.Measurements.insert().values(measurement_type_id=measurement_type_id,
                                                measurement_name=measurement_name,
                                                measurement_description=measurement_description,
                                                measurement_datetime=measurement_datetime)
        ins.compile().params
        result = self.connection.execute(ins)
        self.log_record('Measurement: %s added' % (measurement_name), 'Information')
        return result.inserted_primary_key[0]

    def delete_measurement_datapoints(self, measurement_id):
        delete_query = self.DataPoints.delete().where(self.DataPoints.c.measurement_id == measurement_id)
        self.connection.execute(delete_query)

    def delete_measurement_with_datapoints(self, measurement_id):
        self.delete_measurement_datapoints(measurement_id)
        delete_query = self.Measurements.delete().where(self.Measurements.c.measurement_id == measurement_id)
        self.connection.execute(delete_query)

    def add_datapoint_array(self, point_order, measurement_id, tool_id,
                            point_name_long, point_name_short, point_unit_name,
                            point_value, point_measured):
        four_int_val = number_to_four_ints(point_value)
        self.engine.execute(
            self.DataPoints.insert(),
            [{'point_order': point_order[i],
              'measurement_id': measurement_id,
              'tool_id': tool_id,
              'point_name_long': point_name_long,
              'point_name_short': point_name_short,
              'point_unit_name': point_unit_name,
              'point_sign': four_int_val[i][0],
              'point_mantissa': four_int_val[i][1],
              'point_exponent': four_int_val[i][2],
              'point_bytecount': four_int_val[i][3],
              'point_measured': point_measured,
              } for i in range(point_value.size)])

    def add_datapoint(self, point_order, measurement_id, tool_id,
                      point_name_long, point_name_short, point_unit_name,
                      point_value, point_measured):
        four_int_val = number_to_four_ints(point_value)
        # print 'POINT_VALUE:', point_value
        # print '4 INTS:', four_int_val

        ins = self.DataPoints.insert().values(point_order=point_order, measurement_id=measurement_id, tool_id=tool_id,
                                              point_name_long=point_name_long, point_name_short=point_name_short,
                                              point_unit_name=point_unit_name,
                                              point_sign=four_int_val[0], point_mantissa=four_int_val[1],
                                              point_exponent=four_int_val[2], point_bytecount=four_int_val[3],
                                              point_measured=point_measured)
        ins.compile().params
        result = self.connection.execute(ins)
        return result.inserted_primary_key

    def get_data_points_by_names(self, measurement_id, points_names):
        # t0 = time.time()
        data = {point_name: [] for point_name in points_names}
        # print data
        sel = select([self.DataPoints]).where(self.DataPoints.c.measurement_id == measurement_id).order_by(
            'point_name_short')
        result = self.connection.execute(sel)
        # t1 = time.time()
        # print 'SELECT time', t1-t0, 's'
        for row in result:
            if row['point_name_long'] in points_names:
                sign, mantissa, exponent, bytecount = row['point_sign'], row['point_mantissa'], row['point_exponent'], \
                                                      row['point_bytecount']
                val = np.float(mpmath.mpf((sign, mantissa, exponent, bytecount)))
                data[row['point_name_long']].append(
                    [val, row['point_order'], row['point_measured'], row['point_unit_name']])
        # print data
        for point_name in data.keys():
            data[point_name] = np.array(data[point_name])
            data[point_name] = data[point_name][data[point_name][:, 1].argsort()]
        # t2 = time.time()
        # print 'DATA CONVERSION time', t2-t1, 's'
        # print 'TOTAL time', t2-t0, 's\n'
        return data

    def get_data_points(self, measurement_id, point_name):
        # t0 = time.time()
        sel = select([self.DataPoints]).where(and_(self.DataPoints.c.measurement_id == measurement_id,
                                                   or_(self.DataPoints.c.point_name_long.like(point_name),
                                                       self.DataPoints.c.point_name_short.like(point_name)))).order_by(
            'point_measured')
        result = self.connection.execute(sel)
        # t1 = time.time()
        # print 'SELECT time', t1-t0, 's'
        data_points = []
        measurement_date = []
        point_unit_name = []
        point_order = []
        for row in result:
            sign, mantissa, exponent, bytecount = row['point_sign'], row['point_mantissa'], row['point_exponent'], row[
                'point_bytecount']
            data_points.append(np.float(mpmath.mpf((sign, mantissa, exponent, bytecount))))
            measurement_date.append(row['point_measured'])
            point_unit_name.append(row['point_unit_name'])
            point_order.append(row['point_order'])
        # t2 = time.time()
        # print 'DATA CONVERSION time', t2-t1, 's'
        # print 'TOTAL time', t2-t0, 's\n'
        return np.array(data_points), np.array(point_order), np.array(measurement_date), np.array(point_unit_name)

    def get_next_data_point_order(self, measurement_id, point_name):
        sel = select([func.max(self.DataPoints.c.point_order, type_=Integer)],
                     and_(self.DataPoints.c.measurement_id == measurement_id,
                          or_(self.DataPoints.c.point_name_long.like(point_name),
                              self.DataPoints.c.point_name_short.like(point_name))))
        result = self.connection.execute(sel)
        next_order = result.fetchone()[0]
        next_order = 0 if next_order is None else next_order + 1
        return next_order

    def get_data_point_at_datetime(self, measurement_id, point_name, date_time, interpolation_type='step'):
        if interpolation_type == 'step':
            sel = select([self.DataPoints]).where(and_(and_(self.DataPoints.c.measurement_id == measurement_id,
                                                            or_(self.DataPoints.c.point_name_long.like(point_name),
                                                                self.DataPoints.c.point_name_short.like(point_name))),
                                                       self.DataPoints.c.point_measured <= date_time)).order_by(
                desc(self.DataPoints.c.point_measured))
            result = self.connection.execute(sel)
            row = result.fetchone()
            sign, mantissa, exponent, bytecount = row['point_sign'], row['point_mantissa'], row['point_exponent'], row[
                'point_bytecount']
            data_point = mpmath.mpf((sign, mantissa, exponent, bytecount))
        else:
            data_points, data_points_measured = self.get_data_points(measurement_id, point_name)
            data_points_measured = data_points_measured.astype('datetime64[us]').astype('d')
            Interp = interp1d(data_points_measured, data_points)
            data_point = Interp(np.datetime64(date_time, 'us').astype('d'))
        return data_point, date_time

    def add_monitored_datapoint(self, measurement_type_id, measurement_name, measurement_description,
                                tool_id, point_name_long, point_name_short, point_unit_name, default, val):
        measurement_id = self.add_measurement(measurement_type_id, measurement_name, measurement_description,
                                              measurement_datetime=None, create_new=False)
        if val is None:
            temp, _, _, _ = self.get_data_points(measurement_id, point_name_long)
            if len(temp) > 0:
                ret_val = temp[-1]
            else:
                measurement_id, ret_val = self.add_monitored_datapoint(measurement_type_id, measurement_name,
                                                                       measurement_description, tool_id,
                                                                       point_name_long, point_name_short,
                                                                       point_unit_name, default, default)
        else:
            point_order = self.get_next_data_point_order(measurement_id, point_name_long)
            ret_val = val
            self.add_datapoint(point_order, measurement_id, tool_id,
                               point_name_long, point_name_short, point_unit_name,
                               val, datetime.datetime.now())
        return measurement_id, ret_val
