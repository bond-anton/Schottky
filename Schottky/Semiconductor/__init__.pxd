from BDFunction1D cimport Function
from Schottky.Trap cimport Trap
from Schottky.Helpers.Cache cimport Cache


cdef class Semiconductor(object):
    cdef:
        str __label
        dict __reference
        list __dopants
        Cache __band_gap_cache
        Cache __trap_eq_occupation_cache
        Cache __bulk_charge_cache
        Cache __el_chem_pot_cache
        Function __dopants_n

    cdef double __band_gap_t(self, double temperature)
    cpdef double band_gap_t(self, double temperature)
    cpdef double band_gap_ev_t(self, double temperature)
    cpdef double band_gap_boltzmann_t(self, double temperature)
    cpdef double[:] band_gap(self, double[:] temperature)
    cpdef double[:] band_gap_ev(self, double[:] temperature)
    cpdef double[:] band_gap_boltzmann(self, double[:] temperature)

    cpdef double n_c_t(self, double temperature)
    cpdef double[:] n_c(self, double[:] temperature)
    cpdef double n_v_t(self, double temperature)
    cpdef double[:] n_v(self, double[:] temperature)

    cpdef double pn_t(self, double temperature)
    cpdef double[:] pn(self, double[:] temperature)
    cpdef double n_i_t(self, double temperature)
    cpdef double[:] n_i(self, double[:] temperature)

    cpdef double one_by_n_i_t(self, double temperature)
    cpdef double[:] one_by_n_i(self, double[:] temperature)

    cpdef double n_c_n_i_t(self, double temperature)
    cpdef double[:] n_c_n_i(self, double[:] temperature)
    cpdef double n_v_n_i_t(self, double temperature)
    cpdef double[:] n_v_n_i(self, double[:] temperature)

    cpdef double e_i_t(self, double temperature)
    cpdef double[:] e_i(self, double[:] temperature)
    cpdef double e_i_ev_t(self, double temperature)
    cpdef double[:] e_i_ev(self, double[:] temperature)
    cpdef double e_i_boltzmann_t(self, double temperature)
    cpdef double[:] e_i_boltzmann(self, double[:] temperature)

    cpdef double n_e_t(self, double mu, double temperature)
    cpdef double[:] n_e(self, double[:] mu, double[:] temperature)
    cpdef double n_e_ev_t(self, double mu_ev, double temperature)
    cpdef double[:] n_e_ev(self, double[:] mu_ev, double[:] temperature)
    cpdef double n_e_boltzmann_t(self, double mu, double temperature)
    cpdef double[:] n_e_boltzmann(self, double[:] mu, double[:] temperature)

    cpdef double n_e_to_mu_t(self, double n_e, double temperature)
    cpdef double n_e_to_mu_ev_t(self, double n_e, double temperature)
    cpdef double n_e_to_mu_boltzmann_t(self, double n_e, double temperature)

    cpdef double n_h_t(self, double mu, double temperature)
    cpdef double[:] n_h(self, double[:] mu, double[:] temperature)
    cpdef double n_h_ev_t(self, double mu_ev, double temperature)
    cpdef double[:] n_h_ev(self, double[:] mu_ev, double[:] temperature)
    cpdef double n_h_boltzmann_t(self, double mu, double temperature)
    cpdef double[:] n_h_boltzmann(self, double[:] mu, double[:] temperature)

    cpdef double n_h_to_mu_t(self, double n_h, double temperature)
    cpdef double n_h_to_mu_ev_t(self, double n_h, double temperature)
    cpdef double n_h_to_mu_boltzmann_t(self, double n_h, double temperature)

    cpdef double n_e_n_i_t(self, double mu, double temperature)
    cpdef double[:] n_e_n_i(self, double[:] mu, double[:] temperature)
    cpdef double n_e_n_i_ev_t(self, double mu_ev, double temperature)
    cpdef double[:] n_e_n_i_ev(self, double[:] mu_ev, double[:] temperature)
    cpdef double n_e_n_i_boltzmann_t(self, double mu, double temperature)
    cpdef double[:] n_e_n_i_boltzmann(self, double[:] mu, double[:] temperature)

    cpdef double n_h_n_i_t(self, double mu, double temperature)
    cpdef double[:] n_h_n_i(self, double[:] mu, double[:] temperature)
    cpdef double n_h_n_i_ev_t(self, double mu_ev, double temperature)
    cpdef double[:] n_h_n_i_ev(self, double[:] mu_ev, double[:] temperature)
    cpdef double n_h_n_i_boltzmann_t(self, double mu, double temperature)
    cpdef double[:] n_h_n_i_boltzmann(self, double[:] mu, double[:] temperature)

    cpdef double v_e_t(self, double temperature)
    cpdef double[:] v_e(self, double[:] temperature)
    cpdef double v_h_t(self, double temperature)
    cpdef double[:] v_h(self, double[:] temperature)
    cpdef double trap_eq_occupation(self, Trap trap, double mu, double temperature,
                                    double f_threshold=*, int max_iter=*, bint verbose=*)
    cpdef double bulk_charge(self, double mu, double temperature, double z=*,
                             double f_threshold=*, int max_iter=*, bint verbose=*)
    cpdef double el_chem_pot_t(self, double temperature,
                               double f_threshold=*, int max_iter=*, bint verbose=*)
    cpdef double[:] el_chem_pot(self, double[:] temperature,
                                double f_threshold=*, int max_iter=*, bint verbose=*)
    cpdef double el_chem_pot_ev_t(self, double temperature,
                                  double f_threshold=*, int max_iter=*, bint verbose=*)
    cpdef double[:] el_chem_pot_ev(self, double[:] temperature,
                                   double f_threshold=*, int max_iter=*, bint verbose=*)
    cpdef double work_function_t(self, double temperature,
                                 double f_threshold=*, int max_iter=*, bint verbose=*)
    cpdef double[:] work_function(self, double[:] temperature,
                                  double f_threshold=*, int max_iter=*, bint verbose=*)
    cpdef double work_function_ev_t(self, double temperature,
                                    double f_threshold=*, int max_iter=*, bint verbose=*)
    cpdef double[:] work_function_ev(self, double[:] temperature,
                                     double f_threshold=*, int max_iter=*, bint verbose=*)
    cpdef double work_function_boltzmann_t(self, double temperature,
                                           double f_threshold=*, int max_iter=*, bint verbose=*)
    cpdef double[:] work_function_boltzmann(self, double[:] temperature,
                                            double f_threshold=*, int max_iter=*, bint verbose=*)

    cpdef double mobility_e_point_t(self, double dopants_n, double field, double pn, double temperature)
    cpdef double[:] mobility_e_point(self, double dopants_n, double field, double pn, double[:] temperature)
    cpdef double[:] mobility_e_t(self, double[:] dopants_n, double[:] field, double[:] pn, double temperature)
    cpdef double mobility_h_point_t(self, double dopants_n, double field, double pn, double temperature)
    cpdef double[:] mobility_h_point(self, double dopants_n, double field, double pn, double[:] temperature)
    cpdef double[:] mobility_h_t(self, double[:] dopants_n, double[:] field, double[:] pn, double temperature)

    cpdef double diffusivity_e_point_t(self, double dopants_n, double field, double pn, double temperature)
    cpdef double[:] diffusivity_e_point(self, double dopants_n, double field, double pn, double[:] temperature)
    cpdef double[:] diffusivity_e_t(self, double[:] dopants_n, double[:] field, double[:] pn, double temperature)
    cpdef double diffusivity_h_point_t(self, double dopants_n, double field, double pn, double temperature)
    cpdef double[:] diffusivity_h_point(self, double dopants_n, double field, double pn, double[:] temperature)
    cpdef double[:] diffusivity_h_t(self, double[:] dopants_n, double[:] field, double[:] pn, double temperature)

    cpdef double debye_length_point_t(self, double dopants_n, double temperature)
    cpdef double[:] debye_length_point(self, double dopants_n, double[:] temperature)
    cpdef double[:] debye_length_t(self, double[:] dopants_n, double temperature)
