import numpy as np
from planck import planck_function
from astropy import units as u
from astropy import constants as const

def bb_flux(wavelength, temperature, angular_radius):
   
    f_lambda = (np.pi * u.sr) * planck_function(wavelength, temperature) \
               * (angular_radius)**2

    return f_lambda
