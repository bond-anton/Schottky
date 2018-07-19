from scipy.constants import c, N_A, k, e, m_e, epsilon_0

cdef class Constants(object):

    def __init__(self):
        self.__c = c
        self.__avogadro = N_A
        self.__q = e
        self.__m_e = m_e
        self.__A_R = 1.20173e6  # A*m-2*K-2
        self.__k = k
        self.__epsilon_0 = epsilon_0

    @property
    def c(self):
        return self.__c

    @property
    def avogadro(self):
        return self.__avogadro

    @property
    def q(self):
        return self.__q

    @property
    def m_e(self):
        return self.__m_e

    @property
    def A_R(self):
        return self.__A_R

    @property
    def k(self):
        return self.__k

    @property
    def epsilon_0(self):
        return self.__epsilon_0
