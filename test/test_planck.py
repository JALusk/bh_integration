import unittest
import math
import bh_integration.constants as constants
from bh_integration.planck import planck_function

class TestPlanckFunctionExtrema(unittest.TestCase):

    def setUp(self):
        self.temperature = 5000.0 #Temperature in Kelvin
        self.C1 = 2.0 * constants.h * constants.c**2
        self.C2 = constants.h * constants.c / constants.k

    def test_rayleigh_jeans_law_long_wavelength(self):
        # The Planck function should be within 1% of the Rayleigh-Jeans
        # approximation for lambda * T > 7.7E9 angstrom * K
        wavelength = 7.8E9 / self.temperature
        
        # Expected value is the Rayleigh-Jeans approximation
        expected = self.C1/self.C2 * self.temperature / (wavelength)**(4)
        result = planck_function(wavelength, self.temperature)
        
        self.assertAlmostEqual(expected, result, delta = 0.01 * expected)

    def test_wien_approx_short_wavelength(self):
        # The Planck function should be within 1% of the Wein approximation
        # for lambda * T < 3.0E7 angstrom * K
        wavelength = 2.9E7 / self.temperature
        
        # Expected value is the Wein approximation
        expected = self.C1 * (wavelength)**(-5) * math.exp(-self.C2 / 
                             (wavelength * self.temperature))
        result = planck_function(wavelength, self.temperature)
        
        self.assertAlmostEqual(expected, result, delta = 0.01 * expected)

    def test_rayleigh_jeans_law_short_wavelength(self):
        # The Rayleigh-Jeans law should NOT work for short wavelength
        wavelength = 2.9E7 / self.temperature

        expected = self.C1/self.C2 * self.temperature / (wavelength)**(4)
        result = planck_function(wavelength, self.temperature)

        self.assertNotAlmostEqual(expected, result, delta = 0.01 * expected)

    def test_wien_approx_long_wavelength(self):
        # The Wien approximation should NOT work for long wavelength
        wavelength = 7.8E9 / self.temperature

        expected = self.C1 * (wavelength)**(-5) * math.exp(-self.C2 / 
                             (wavelength * self.temperature))
        result = planck_function(wavelength, self.temperature)

        self.assertNotAlmostEqual(expected, result, delta = 0.01 * expected)
