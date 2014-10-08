import bh_integration.constants as constants
import math

def planck_function(wavelength, temperature):

    a = 2.0 * constants.h * constants.c**2
    b = constants.h * constants.c / constants.k

    B_lambda = (a / wavelength**5)/(math.exp(b / (wavelength * temperature)) - 1)
    
    return B_lambda
