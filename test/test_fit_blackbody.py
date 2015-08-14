import unittest
import numpy as np
from astropy import units as u
from scipy.optimize import curve_fit
from bh_integration.planck import planck_function
from bh_integration.fit_blackbody import *

class TestFitBlackbody(unittest.TestCase):

    def setUp(self):
        self.wavelength = 5000. * u.Angstrom
        self.temperature = 13000. * u.K
        self.angular_radius = 0.2e-10
        self.flux_array = u.Quantity([5.87565760e-12, 3.79417256e-11, 
                                      6.40953095e-11, 5.90280490e-11, 
                                      3.41930932e-11], 
                                      u.erg / (u.Angstrom * u.cm**2 * u.s))
        self.flux_uncertainties = u.Quantity([1.0e-13, 1.0e-12, 
                                      1.0e-12, 1.0e-12, 
                                      1.0e-12], 
                                      u.erg / (u.Angstrom * u.cm**2 * u.s))

        self.eff_wl_array = u.Quantity([3660., 4380., 5450., 6410., 7980.],
                                        u.Angstrom)

    def test_bb_flux_returns_expected_flux(self):
        expected = (np.pi * u.sr) \
                   * planck_function(self.wavelength, self.temperature) \
                   * (self.angular_radius)**2
        result = bb_flux(self.wavelength, self.temperature, 
                         self.angular_radius)

        self.assertEqual(expected, result)
    
    def test_bb_flux_nounits_returns_expected_flux(self):
        expected = (np.pi * u.sr) \
                   * planck_function(self.wavelength, self.temperature) \
                   * (self.angular_radius)**2
        result = bb_flux_nounits(self.wavelength, self.temperature, 
                         self.angular_radius)

        self.assertEqual(expected.value, result)

    def test_bb_fit_parameters_returns_expected_parameters(self):
        popt, pcov = curve_fit(bb_flux_nounits, self.eff_wl_array, self.flux_array,
                               p0=[5000, 1.0e-10])
        expected_temp = u.Quantity(popt[0], unit=u.K)
        expected_radius = popt[1]
        expected_chisq = calculate_chisq(self.flux_array, self.flux_uncertainties, self.eff_wl_array, bb_flux, popt)
        result_temp, result_radius, result_chisq = bb_fit_parameters(self.eff_wl_array,
                                                       self.flux_array, self.flux_uncertainties)
        self.assertEqual((expected_temp, expected_radius, expected_chisq), 
                         (result_temp, result_radius, result_chisq))
class TestChiSquared(unittest.TestCase):

    def setUp(self):
        self.y_data = np.array([1, 2, 3, 4.5, 5])
        self.y_data_uncertainty = np.array([0.1, 0.1, 0.1, 0.1, 0.1])
        self.x_data = np.array([1, 2, 3, 4, 5])
        self.parameters = [1,0]
        
    def test_chi_squared_calculation(self):
        expected = 25.0

        def func(x, *parameters):
            m = parameters[0]
            b = parameters[1]
            return m * x + b

        result = calculate_chisq(self.y_data, self.y_data_uncertainty, self.x_data, func, self.parameters)

        self.assertEqual(expected, result)

    def test_chi_squared_quadratic(self):
        expected = 57225.0

        def func(x, *parameters):
            a = parameters[0]
            b = 0
            c = parameters[1]

            return a * x**2 + b * x + c

        result = calculate_chisq(self.y_data, self.y_data_uncertainty, self.x_data, func, self.parameters)
        self.assertEqual(expected, result)
