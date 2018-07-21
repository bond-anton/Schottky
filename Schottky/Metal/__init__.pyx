from __future__ import division, print_function

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
        return self.__work_function / constant.__q

    @work_function_ev.setter
    def work_function_ev(self, double work_function_ev):
        self.__work_function = work_function_ev * constant.__q

    def __str__(self):
        return 'Metal: %s, Workfunction: %2.2f eV (%2.2g J)' % (self.label,
                                                                self.work_function_ev,
                                                                self.work_function)
