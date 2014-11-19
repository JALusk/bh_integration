import numpy as np
from astropy import units as u
import scipy.integrate as integrate
from bh_integration.mag2flux import *
from bh_integration.fqbol import integrate_fqbol
from bh_integration.fit_blackbody import *

def ccm89_deredden_magnitudes(key, magnitudes, av):
    mags_dereddened = np.array([])
    
    for i, filter_band in enumerate(key):
        filter_mean_extinction = get_filter_mean_extinction(filter_band)
        dereddened_mag = magnitudes[i] - av * filter_mean_extinction
        mags_dereddened = np.append(mags_dereddened, dereddened_mag)

    return mags_dereddened

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
    magnitudes = ccm89_deredden_magnitudes(key, magnitudes, av)
    
    wavelength_array, flux_array = build_flux_wl_array(key, magnitudes)
    
    shortest_wl = np.amin(wavelength_array)
    longest_wl = np.amax(wavelength_array)

    fqbol = integrate_fqbol(wavelength_array, flux_array)
    
    temperature, angular_radius = bb_fit_parameters(wavelength_array, 
                                                    flux_array)
    ir_correction = calculate_ir_correction(bb_flux(wavelength, temperature,
                                                    angular_radius),
                                                    longest_wl)
    uv_correction = calculate_uv_correction_using_blackbody(bb_flux(wavelength, 
                                                            temperature, 
                                                            angular_radius), 
                                                            shortest_wl)
    fbol = fqbol + ir_correction + uv_correction
    return fbol
