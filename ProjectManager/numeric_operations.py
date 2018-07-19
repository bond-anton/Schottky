from __future__ import division

import mpmath as mp

import numpy as np

def number_to_four_ints(x):
    '''
    returns array of four ints
    sign, mantissa, exponent, bytecount
    '''
    if isinstance(x, np.ndarray):
        return np.array([number_to_four_ints(x_i) for x_i in x])
    else:
        g = mp.mpf(x)
        sign, mantissa, exponent, bc = g._mpf_
        a = np.array([int(sign), int(mantissa), int(exponent), int(bc)], dtype=np.int)
        print(a)
        return a