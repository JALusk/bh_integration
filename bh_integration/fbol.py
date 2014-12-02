import numpy as np
from astropy import units as u
import scipy.integrate as integrate
from mag2flux import *
from fqbol import integrate_fqbol
from specutils import extinction
from fit_blackbody import *

def build_flux_wl_array(key, magnitudes):
    fluxes = np.array([])
    wavelengths = np.array([])

    for i, filter_band in enumerate(key):
        flux, effective_wl = mag2flux(filter_band, magnitudes[i])
        fluxes = np.append(fluxes, flux)
        wavelengths = np.append(wavelengths, effective_wl)

    return u.Quantity(wavelengths, effective_wl.unit), u.Quantity(fluxes, flux.unit)

def integrate_fqbol(wavelength_array, flux_array):
    return np.trapz(flux_array, wavelength_array)

def ir_correction(temperature, angular_radius, longest_wl):
    ir_correction = integrate.quad(bb_flux_nounits, longest_wl, np.inf,
                                   args=(temperature, angular_radius)) 
    return ir_correction

def uv_correction_blackbody(temperature, angular_radius, shortest_wl):
    uv_correction = integrate.quad(bb_flux_nounits, 0, shortest_wl, 
                                   args=(temperature, angular_radius))
    return uv_correction

def uv_correction_linear(shortest_wl, shortest_flux):
    fluxes = u.Quantity([0.0, shortest_wl.value], shortest_wl.unit)
    wavelengths = u.Quantity([0.0, shortest_flux.value], shortest_flux.unit)
    uv_correction = np.trapz(fluxes, wavelengths)
    return uv_correction

def calculate_fbol(key, magnitudes, av):
    wavelength_array, flux_array = build_flux_wl_array(key, magnitudes)

    flux_array_unred = flux_array * extinction.reddening(wavelength_array, av, model='od94')
    
    shortest_wl = np.amin(wavelength_array)
    longest_wl = np.amax(wavelength_array)

    fqbol = integrate_fqbol(wavelength_array, flux_array)
    
    temperature, angular_radius = bb_fit_parameters(wavelength_array, 
                                                    flux_array)
    ir_corr = ir_correction(temperature.value, angular_radius, longest_wl.value)[0]
    uv_corr = uv_correction_blackbody(temperature.value, angular_radius, 
                                            shortest_wl.value)[0]
    fbol = fqbol.value + ir_corr + uv_corr
    return fbol
