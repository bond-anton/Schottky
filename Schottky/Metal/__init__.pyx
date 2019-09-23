from cython cimport boundscheck, wraparound
from cpython.array cimport array, clone

from Schottky.Constants cimport constant


cdef class Metal(object):
    '''
    Metal electrode class
    '''

    def __init__(self, str label, double work_function):
        '''
        Constructor
        '''
        self.__label = label
        self.__work_function = work_function

    @property
    def label(self):
        return self.__label

    @label.setter
    def label(self, str label):
        self.__label = label

    @property
    def work_function(self):
        return self.__work_function

    @work_function.setter
    def work_function(self, double work_function):
        self.__work_function = work_function

    @property
    def work_function_ev(self):
        return constant.joule_to_ev_point(self.__work_function)

    @work_function_ev.setter
    def work_function_ev(self, double work_function_ev):
        self.__work_function = constant.ev_to_joule_point(work_function_ev)

    cpdef double work_function_boltzmann_t(self, double temperature):
        return constant.joule_to_boltzmann_point(self.__work_function, temperature)

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] work_function_boltzmann(self, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.work_function_boltzmann_t(temperature[i])
        return result

    def __str__(self):
        return 'Metal: %s, Workfunction: %2.2f eV (%2.2g J)' % (self.label,
                                                                self.work_function_ev,
                                                                self.work_function)
