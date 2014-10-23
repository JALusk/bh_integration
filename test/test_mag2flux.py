import unittest
import numpy as np
from astropy import constants as const
from astropy import units as u
from bh_integration.mag2flux import *
from yaml import load

stream = file('filter_data.yaml', 'r')
filter_data = load(stream)

class TestMagnitudeToFluxConversion(unittest.TestCase):

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

    def test_flux_at_mag_zero(self):
        mag = 0.0
        expected = self.flux_at_zero_mag
        result = mag2flux(self.filter_band, 0.0)
        
        self.assertEqual(expected.value, result.value)

    def test_get_filter_parameters_bad_filter_name(self):
        self.assertRaises(KeyError, get_filter_parameters, 'u')

    def test_get_filter_parameters_bad_argument_int(self):
        self.assertRaises(TypeError, get_filter_parameters, 12)
