from BDSpace.Field.Field cimport Field
from BDSpace.Field.SuperposedField cimport SuperposedField
from Schottky.Potential.ExternalField cimport ExternalField
from Schottky.Trap cimport Trap


cdef class TrapPotential(SuperposedField):
    cdef:
        Trap __trap
        Field __trap_field
        ExternalField __external_field
    cpdef double emission_rate_enhancement(self, double temperature=*)


cdef class NullPotential(TrapPotential):
    pass


cdef class PointLikeInExternalField(TrapPotential):
    cdef:
        double __r_min
        double __r_max
        int __phi_resolution
        int __theta_resolution
    cpdef double max_energy_r_point(self, double theta, double phi)
    cpdef double[:] max_energy_r(self, double theta, double[:] phi)
    cpdef double energy_lowering_point(self, double theta, double phi)
    cpdef double[:] energy_lowering_phi_range(self, double theta, double[:] phi)
    cpdef double[:] energy_lowering_theta_range(self, double phi, double[:] theta)


cdef class SphericallySymmetricInExternalField(PointLikeInExternalField):
    cpdef bint is_aligned(self)


cdef class HyperbolicInExternalField(SphericallySymmetricInExternalField):
    pass
