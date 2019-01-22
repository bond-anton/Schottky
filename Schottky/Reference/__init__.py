from Schottky import constant


q = constant.q

database = [
    {
        'label': 'Si',
        'Z': 14,  # mean atomic number
        'density': 2329.0,  # kg/m^3
        'epsilon': 11.8,  # 'mu' : 0.1,
        'thermionic_emission': {
            'A_R_coeff_n': 2.1,
            'A_R_coeff_p': 0.6
        },
        'N_c0': 6.2e21,
        'N_v0': 3.5e21,
        'affinity': q * 4.05,
        'band_gap': {
            'E_g0': q * 1.17,
            'alpha': q * 4.73e-4,
            'beta': 636.0
        },
        'carrier_mass': {
            'm_e_coeff': 0.19,
            'm_h_coeff': 0.16
        },
        'mobility_e': {
            'mu_L0': 1430,  # cm2/(V*s)
            'alpha': 2.2,
            'A': 4.61e17,  # cm-1 * V-1 * s-1 * K-3/2
            'B': 1.52e15,  # cm-3 * K-2
            'v_s': 1.0e7,  # cm/s
            'beta': 2,
        },
        'mobility_h': {
            'mu_L0': 495,  # cm2/(V*s)
            'alpha': 2.2,
            'A': 1.0e17,  # cm-1 * V-1 * s-1 * K-3/2
            'B': 6.25e14,  # cm-3 * K-2
            'v_s': 1.0e7,  # cm/s
            'beta': 1,
        },
    },
]
