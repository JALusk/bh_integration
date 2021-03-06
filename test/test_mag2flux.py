import unittest
import numpy as np
from astropy import constants as const
from astropy import units as u
from bh_integration.mag2flux import *
from yaml import load

stream = file('filter_data.yaml', 'r')
filter_data = load(stream)

class TestGetFilterParameters(unittest.TestCase):

    def setUp(self):
        self.filter_band = "V"
        self.effective_wl = filter_data[self.filter_band]['effective_wl'] * u.AA
        self.flux_at_zero_mag = filter_data[self.filter_band]['flux_at_zero_mag'] * \
                                (u.erg / (u.s * u.cm**2 * u.AA))

    def test_get_filter_parameters_returns_effective_wl(self):
        expected = self.effective_wl
        result = get_filter_parameters(self.filter_band)[0]

        self.assertEqual(expected.value, result.value)

    def test_get_filter_parameters_returns_flux_at_zero_mag(self):
        expected = self.flux_at_zero_mag
        result = get_filter_parameters(self.filter_band)[1]

        self.assertEqual(expected.value, result.value)

   
    def test_get_filter_parameters_bad_filter_name(self):
        self.assertRaises(KeyError, get_filter_parameters, 'u')

    def test_get_filter_parameters_bad_argument_int(self):
        self.assertRaises(TypeError, get_filter_parameters, 12)
        
class TestMag2Flux(unittest.TestCase):

    def setUp(self):
        self.filter_band = "V"
        self.magnitude = 8.8
        self.uncertainty = 0.02
        self.effective_wl = filter_data[self.filter_band]['effective_wl'] * u.AA
        self.flux_at_zero_mag = filter_data[self.filter_band]['flux_at_zero_mag'] * \
                                (u.erg / (u.s * u.cm**2 * u.AA))
    
    def test_mag2flux_converts_mag_to_correct_flux(self):
        expected = self.flux_at_zero_mag * 10**(-0.4 * self.magnitude)
        result_flux, result_uncertainty, result_wl = mag2flux(self.filter_band, 
                                                              self.magnitude, 
                                                              self.uncertainty)

        self.assertEqual(expected.value, result_flux.value)

    def test_mag2flux_converts_mag_to_correct_flux_uncertainty(self):
        expected = np.sqrt((self.flux_at_zero_mag * -0.4 * np.log(10) * 10**(-0.4 * self.magnitude) * self.uncertainty)**2)
        result_flux, result_uncertainty, result_wl = mag2flux(self.filter_band, 
                                                              self.magnitude, 
                                                              self.uncertainty)

        self.assertAlmostEqual(expected.value, result_uncertainty.value)


    def test_mag2flux_returns_correct_effective_wl(self):
        expected = self.effective_wl
        result_flux, result_uncertainty, result_wl = mag2flux(self.filter_band, self.magnitude,
                                          self.uncertainty)

        self.assertEqual(expected.value, result_wl.value)

    def test_flux_at_mag_zero(self):
        mag = 0.0
        expected = self.flux_at_zero_mag
        result_flux, result_uncertainty, result_wl = mag2flux(self.filter_band, 0.0, self.uncertainty)
        
        self.assertEqual(expected.value, result_flux.value)
