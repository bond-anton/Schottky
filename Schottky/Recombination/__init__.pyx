from BDFunction1D.Standard cimport Constant, Zero


cdef class ConstantRecombinationFunction(Constant):

    def __init__(self, double rate):
        super(ConstantRecombinationFunction, self).__init__(rate)

    @property
    def rate(self):
        return self.c

    @rate.setter
    def rate(self, double rate):
        self.c = rate


cdef class ZeroRecombinationFunction(Zero):

    @property
    def rate(self):
        return 0.0
