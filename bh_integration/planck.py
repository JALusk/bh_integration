import bh_integration.constants as constants
import math

def planck_function(wavelength, temperature):

    C1 = 2.0 * constants.h * constants.c**2
    C2 = constants.h * constants.c / constants.k

    B_lambda = (C1 / wavelength**5)/(math.exp(C2 / (wavelength * temperature)) - 1)
    
    return B_lambda
