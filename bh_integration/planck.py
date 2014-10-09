import numpy as np
from astropy import constants as const
from astropy import units as u

def planck_function(wavelength, temperature):

    C1 = 2.0 * const.h.cgs * const.c.cgs**2
    C2 = const.h.cgs * const.c.cgs / const.k_B.cgs

    B_lambda = (C1 / wavelength**5)/(np.exp(C2 / (wavelength * temperature)) - 1) / u.sr
    B_lambda = B_lambda.to(u.erg / (u.s * u.cm**2 * u.AA * u.sr))

    return B_lambda
