import numpy as np
from astropy import units as u
from yaml import load

stream = file('filter_data.yaml', 'r')
filter_data = load(stream)

def get_filter_parameters(filter_band):
    """Fetches effective wavelength and flux zeropoint of a filter.

    Args:
        filter_band: StringType - Name of the filter band (ex: 'U').

    Returns:
        A tuple containing two astropy quantities: the effective wavelength
        in Angstroms, and the flux at zero magnitudes in erg/s/cm^s/A.

        (effective_wl, flux_at_zero_mag).

    Raises:
        TypeError: The argument given is not a string.
        KeyError: The argument given is not a valid filter name.
    """
    if not isinstance(filter_band, str):
        raise TypeError('filter_band must be StringType')
    elif filter_band in filter_data:
        effective_wl = filter_data[filter_band]['effective_wl']
        flux_at_zero_mag = filter_data[filter_band]['flux_at_zero_mag']
    else:
        raise KeyError('Your specified filter (%s) is not available', filter_band)
   
    # Attach CGS units to parameters
    effective_wl = effective_wl * u.AA
    flux_at_zero_mag = flux_at_zero_mag * (u.erg / (u.s * u.cm**2 * u.AA)) 

    return effective_wl, flux_at_zero_mag

def mag2flux(filter_band, magnitude):
    """Converts an observed magnitude in a filter band to an average flux.

    Args:
        filter_band: StringType - Name of the filter band (ex: 'U').
        magnitude: FloatTypea - Apparent magnitude in the filter band.

    Returns:
        a tuple of two astropy quantities, the flux in erg/s/cm^2/A. and
        the effective wavelength of the filter in Angstrom

        (flux, effective_wl)
    """
    effective_wl, flux_at_zero_mag = get_filter_parameters(filter_band)
    flux = flux_at_zero_mag * 10**(-0.4 * magnitude)

    return flux, effective_wl
