import numpy as np
from scipy.optimize import curve_fit
from astropy import units as u
from planck import planck_function

def bb_flux(wavelength, temperature, angular_radius):
    bb_flux = (np.pi * u.sr) * planck_function(wavelength, temperature) * (angular_radius)**2

    return bb_flux

def bb_flux_nounits(wavelength, temperature, angular_radius):
    flux = bb_flux(wavelength, temperature, angular_radius)

    return flux.value

    
def bb_fit_parameters(wavelength_array, flux_array):
    popt, pcov = curve_fit(bb_flux, wavelength_array, flux_array, p0=[5000, 1.0e-10])
    temperature = u.Quantity(popt[0], unit=u.K)
    angular_radius = popt[1]

    return temperature, angular_radius
