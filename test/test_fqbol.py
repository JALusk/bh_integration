import unittest
import numpy as np
from astropy import constants as const
from astropy import units as u
from bh_integration.mag2flux import *
from bh_integration.fqbol import *

class TestBuildFluxWLArray(unittest.TestCase):

    def setUp(self):
        # Test data from SN 1987A, JD 2446920.495
        # WARNING: NOT CORRECTED FOR REDDENING
        self.key = ['U', 'B']
        self.magnitudes = [7.192, 4.641]

    def test_build_flux_wl_array_correct_flux_U(self):
        expected_flux, expected_wl = mag2flux(self.key[0], self.magnitudes[0])
        fluxes, wavelengths = build_flux_wl_array(self.key, self.magnitudes)
        result_flux = fluxes[0]
        self.assertEqual(expected_flux.value, result_flux.value)

    def test_build_flux_wl_array_correct_flux_B(self):
        expected_flux, expected_wl = mag2flux(self.key[1], self.magnitudes[1])
        fluxes, wavelengths = build_flux_wl_array(self.key, self.magnitudes)
        result_flux = fluxes[1]
        self.assertEqual(expected_flux.value, result_flux.value)

    def test_build_flux_wl_array_correct_wl_U(self):
        expected_flux, expected_wl = mag2flux(self.key[0], self.magnitudes[0])
        fluxes, wavelengths = build_flux_wl_array(self.key, self.magnitudes)
        result_wl = wavelengths[0]
        self.assertEqual(expected_wl.value, result_wl.value)

    def test_build_flux_wl_array_correct_wl_B(self):
        expected_flux, expected_wl = mag2flux(self.key[1], self.magnitudes[1])
        fluxes, wavelengths = build_flux_wl_array(self.key, self.magnitudes)
        result_wl = wavelengths[1]
        self.assertEqual(expected_wl.value, result_wl.value)

    def test_integrate_fqbol(self):
        U_flux, U_wl = mag2flux(self.key[0], self.magnitudes[0])
        B_flux, B_wl = mag2flux(self.key[1], self.magnitudes[1])
        expected = np.trapz([U_flux.value, B_flux.value], [U_wl.value, B_wl.value])
        flux_array, magnitude_array = build_flux_wl_array(self.key, self.magnitudes)
        result = integrate_fqbol(flux_array, magnitude_array)
        self.assertEqual(expected, result.value)
