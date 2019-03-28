from libc.math cimport M_PI, abs, fabs, sqrt, sin, cos, exp
from cython cimport boundscheck, wraparound
from cpython.array cimport array, clone
from BDSpace.Field.Field cimport Field, ConstantScalarConservativeField
from BDSpace.Field.SuperposedField cimport SuperposedField
from BDSpace.Field.SphericallySymmetric cimport SphericallySymmetric, HyperbolicPotentialSphericalConservativeField
from Schottky.Potential.ExternalField cimport ExternalField
from Schottky.Trap cimport Trap, NullTrap
from Schottky.Constants cimport constant
from ._helpers cimport trapz_1d, linspace


cdef class TrapPotential(SuperposedField):

    def __init__(self, str name, Field trap_field, ExternalField external_field=None, Trap trap=None):
        cdef:
            array[double] direction
        if trap is None:
            self.__trap = NullTrap()
        else:
            self.__trap = trap
        self.__trap_field = trap_field
        if external_field is None:
            direction = clone(array('d'), 3, zero=False)
            direction[0] = 0.0
            direction[1] = 0.0
            direction[2] = 1.0
            self.__external_field = ExternalField(name='External Field',
                                                  direction=direction,
                                                  magnitude=0.0)
        else:
            self.__external_field = external_field
        super(TrapPotential, self).__init__(name, [self.__trap_field, self.__external_field])

    @property
    def trap(self):
        return self.__trap

    @trap.setter
    def trap(self, Trap trap):
        self.__trap = trap

    cpdef double emission_rate_enhancement(self, double temperature=300, double f=0.0):
        return 1.0


cdef class NullPotential(TrapPotential):

    def __init__(self, str name, ExternalField external_field=None, Trap trap=None):
        super(NullPotential, self).__init__(name,
                                            ConstantScalarConservativeField(name='Null Potential',
                                                                            field_type='Electric Field',
                                                                            potential=0.0),
                                            external_field, trap)

    cpdef double emission_rate_enhancement(self, double temperature=300, double f=0.0):
        return 1.0


cdef class PointLikeInExternalField(TrapPotential):

    def __init__(self, str name, Field point_like, ExternalField external_field=None, Trap trap=None,
                 double r_min=1.0e-11, double r_max=1.0e-5,
                 int phi_resolution=10, int theta_resolution=50):
        if fabs(r_min) < fabs(r_max):
            self.__r_min = fabs(r_min)
            self.__r_max = fabs(r_max)
        else:
            self.__r_max = fabs(r_min)
            self.__r_min = fabs(r_max)
        if self.__r_min == self.__r_max:
            if self.__r_min == 0.0:
                self.__r_max = 1.0e-5
            else:
                self.__r_max = 2 * self.__r_min
        if abs(phi_resolution) > 0:
            self.__phi_resolution = abs(phi_resolution)
        else:
            self.__phi_resolution = 10
        if abs(theta_resolution) > 0:
            self.__theta_resolution = abs(theta_resolution)
        else:
            self.__theta_resolution = 50
        super(PointLikeInExternalField, self).__init__(name, point_like, external_field, trap)

    @property
    def r_min(self):
        return self.__r_min

    @r_min.setter
    def r_min(self, double r_min):
        if fabs(r_min) < self.__r_max:
            self.__r_min = fabs(r_min)

    @property
    def r_max(self):
        return self.__r_max

    @r_max.setter
    def r_max(self, double r_max):
        if fabs(r_max) > self.__r_min:
            self.__r_max = fabs(r_max)

    @property
    def phi_resolution(self):
        return self.__phi_resolution

    @phi_resolution.setter
    def phi_resolution(self, int phi_resolution):
        if abs(phi_resolution) > 0:
            self.__phi_resolution = abs(phi_resolution)

    @property
    def theta_resolution(self):
        return self.__theta_resolution

    @theta_resolution.setter
    def theta_resolution(self, int theta_resolution):
        if abs(theta_resolution) > 0:
            self.__theta_resolution = abs(theta_resolution)

    @boundscheck(False)
    @wraparound(False)
    cpdef double max_energy_r_point(self, double theta, double phi):
        cdef:
            double gr = (sqrt(5.0) + 1.0) / 2.0
            double a = self.__r_min, b = self.__r_max
            double tol = 1.6e-19
            array[double] template = array('d')
            array[double] point_c = clone(template, 3, zero=False)
            array[double] point_d = clone(template, 3, zero=False)
            double c_val, d_val
        point_c[0] = b - (b - a) / gr
        point_d[0] = a + (b - a) / gr
        point_c[1] = theta
        point_d[1] = theta
        point_c[2] = phi
        point_d[2] = phi
        while fabs(point_c[0] - point_d[0]) > tol:
            c_val = fabs(self.scalar_field_polar_point(point_c))
            d_val = fabs(self.scalar_field_polar_point(point_d))
            if c_val < d_val:
                b = point_d[0]
            else:
                a = point_c[0]
            # we recompute both c and d here to avoid loss of precision,
            # which may lead to incorrect results or infinite loop
            point_c[0] = b - (b - a) / gr
            point_d[0] = a + (b - a) / gr
        return (b + a) / 2

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] max_energy_r(self, double theta, double[:] phi):
        cdef:
            int s = phi.shape[0]
            array[double] result = clone(array('d'), s, zero=False)
        for i in range(s):
            result[i] = self.max_energy_r_point(theta, phi[i])
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double energy_lowering_point(self, double theta, double phi):
        cdef:
            array[double] rtp = clone(array('d'), 3, zero=False)
        rtp[0] = self.max_energy_r_point(theta, phi)
        rtp[1] = theta
        rtp[2] = phi
        if fabs(rtp[0] - self.__r_max) > 1e-15:
            return fabs(self.scalar_field_polar_point(rtp))
        else:
            return 0.0

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] energy_lowering_phi_range(self, double theta, double[:] phi):
        cdef:
            int s = phi.shape[0]
            array[double] result = clone(array('d'), s, zero=False)
        for i in range(s):
            result[i] = self.energy_lowering_point(theta, phi[i])
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] energy_lowering_theta_range(self, double phi, double[:] theta):
        cdef:
            int s = theta.shape[0]
            array[double] result = clone(array('d'), s, zero=False)
        for i in range(s):
            result[i] = self.energy_lowering_point(theta[i], phi)
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double emission_rate_enhancement(self, double temperature=300, double f=0.0):
        cdef:
            int i, j
            array[double] template = array('d')
            double[:] phi, theta, integrand_theta, integrand_phi
        integrand_phi = clone(template, self.__phi_resolution, zero=False)
        integrand_theta = clone(template, self.__theta_resolution, zero=False)
        phi = linspace(0.0, 2 * M_PI, self.__phi_resolution)
        theta = linspace(0.0, M_PI, self.__theta_resolution)
        for i in range(self.__phi_resolution):
            for j in range(self.__theta_resolution):
                integrand_theta[j] = exp(
                    -self.energy_lowering_point(theta[j], phi[i]) / constant.__k / temperature * constant.__q
                ) * sin(theta[j])
            integrand_phi[i] = trapz_1d(integrand_theta, theta)
        return trapz_1d(integrand_phi, phi)


cdef class SphericallySymmetricInExternalField(PointLikeInExternalField):

    def __init__(self, str name, SphericallySymmetric point_like, ExternalField external_field=None, Trap trap=None,
                 double r_min=1.0e-11, double r_max=1.0e-5,
                 int phi_resolution=10, int theta_resolution=50):
        cdef:
            double[:] rot_axis = clone(array('d'), 3, zero=True)
        super(SphericallySymmetricInExternalField, self).__init__(name, point_like, external_field, trap,
                                                                  r_min, r_max,
                                                                  phi_resolution=phi_resolution,
                                                                  theta_resolution=theta_resolution)

    cpdef bint is_aligned(self):
        cdef:
            bint result
        result = fabs(self.__external_field.direction[0]) < 1.0e-10
        result &= fabs(self.__external_field.direction[1]) < 1.0e-10
        result &= fabs(self.__external_field.direction[2] - 1.0) < 1.0e-10
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] max_energy_r(self, double theta, double[:] phi):
        cdef:
            int s = phi.shape[0]
            array[double] result = clone(array('d'), s, zero=False)
            bint aligned = self.is_aligned()
        result[0] = self.max_energy_r_point(theta, phi[0])
        if aligned:
            for i in range(1, s):
                result[i] = result[0]
        else:
            for i in range(1, s):
                result[i] = self.max_energy_r_point(theta, phi[i])
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] energy_lowering_phi_range(self, double theta, double[:] phi):
        cdef:
            int s = phi.shape[0]
            array[double] result = clone(array('d'), s, zero=False)
            bint aligned = self.is_aligned()
        result[0] = self.energy_lowering_point(theta, phi[0])
        if aligned:
            for i in range(1, s):
                result[i] = result[0]
        else:
            for i in range(1, s):
                result[i] = self.energy_lowering_point(theta, phi[i])
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double emission_rate_enhancement(self, double temperature=300, double f=0.0):
        cdef:
            int i, j
            array[double] template = array('d')
            double[:] phi, theta, integrand_theta, integrand_phi
            bint aligned = self.is_aligned()
        integrand_theta = clone(template, self.__theta_resolution, zero=False)
        theta = linspace(0.0, M_PI, self.__theta_resolution)
        if aligned:
            for j in range(self.__theta_resolution):
                integrand_theta[j] = exp(
                    -self.energy_lowering_point(theta[j], 0.0) / constant.__k / temperature * constant.__q
                ) * sin(theta[j])
            return trapz_1d(integrand_theta, theta) * M_PI * 2
        else:
            integrand_phi = clone(template, self.__phi_resolution, zero=False)
            phi = linspace(0.0, 2 * M_PI, self.__phi_resolution)
            for i in range(self.__phi_resolution):
                for j in range(self.__theta_resolution):
                    integrand_theta[j] = exp(
                        -self.energy_lowering_point(theta[j], phi[i]) / constant.__k / temperature * constant.__q
                    ) * sin(theta[j])
                integrand_phi[i] = trapz_1d(integrand_theta, theta)
            return trapz_1d(integrand_phi, phi)


cdef class HyperbolicInExternalField(SphericallySymmetricInExternalField):

    def __init__(self, str name, HyperbolicPotentialSphericalConservativeField point_charge,
                 ExternalField external_field=None, Trap trap=None,
                 double r_min=1.0e-11, double r_max=1.0e-5,
                 int phi_resolution=10, int theta_resolution=50):
        cdef:
            double[:] rot_axis = clone(array('d'), 3, zero=True)
        super(HyperbolicInExternalField, self).__init__(name, point_charge, external_field, trap,
                                                        r_min, r_max,
                                                        phi_resolution=phi_resolution,
                                                        theta_resolution=theta_resolution)

    @boundscheck(False)
    @wraparound(False)
    cpdef double max_energy_r_point(self, double theta, double phi):
        cdef:
            double r, f = self.__external_field.magnitude
            bint aligned = self.is_aligned()
        if aligned:
            if f > 1.0e-10 and theta < M_PI / 2.0:
                return sqrt(self.__trap_field.a / f / cos(theta))
            else:
                return self.__r_max
        return super(HyperbolicInExternalField, self).max_energy_r_point(theta, phi)
