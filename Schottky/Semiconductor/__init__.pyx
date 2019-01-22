from __future__ import division, print_function

from cython cimport boundscheck, wraparound

from Schottky.Dopant cimport Dopant


cdef class Semiconductor(object):

    def __init__(self, str label, dict reference, list dopants=None):
        self.__label = label
        self.__reference = reference
        self.__dopants = []
        if dopants is not None:
            for dopant in dopants:
                if isinstance(dopant, Dopant):
                    self.__dopants.append(dopant)

    @property
    def label(self):
        return self.__label

    @label.setter
    def label(self, str label):
        self.__label = label

    @property
    def reference(self):
        return self.__reference

    @reference.setter
    def reference(self, dict reference):
        self.__reference = reference

    @property
    def dopants(self):
        return self.__dopants

    @dopants.setter
    def dopants(self, list dopants):
        self.__dopants = []
        for dopant in dopants:
            if isinstance(dopant, Dopant):
                self.__dopants.append(dopant)
