import numpy as np
from astropy import units as u
from bh_integration.mag2flux import *

def build_flux_wl_array(key, magnitudes):
    fluxes = np.array([])
    wavelengths = np.array([])

    for i, filter_band in enumerate(key):
        flux, effective_wl = mag2flux(filter_band, magnitudes[i])
        fluxes = np.append(fluxes, flux)
        wavelengths = np.append(wavelengths, effective_wl)

    return u.Quantity(fluxes, flux.unit), u.Quantity(wavelengths, effective_wl.unit)

def integrate_fqbol(fluxes, wavelengths):
    return np.trapz(fluxes, wavelengths)
