import unittest
import numpy as np
from astropy import units as u
import scipy.integrate as integrate
from bh_integration.mag2flux import *
from bh_integration.fbol import *
from bh_integration.fit_blackbody import *

class TestFbol(unittest.TestCase):
    
    def setUp(self):
        self.key = np.array(['U', 'B', 'V', 'R', 'I'])
        self.mean_extinction = np.array([1.569, 1.337, 1.000, 0.751, 0.479])
        self.magnitudes = np.array([7.129, 5.554, 4.383, 3.917, 3.794])
        self.av = 0.435
        self.flux_array = u.Quantity([5.87565760e-12, 3.79417256e-11, 
                                      6.40953095e-11, 5.90280490e-11, 
                                      3.41930932e-11], 
                                      u.erg / (u.Angstrom * u.cm**2 * u.s))
        self.eff_wl_array = u.Quantity([3660., 4380., 5450., 6410., 7980.],
                                        u.Angstrom)
    
    def test_ccm89_deredden_magnitudes(self):
        expected = np.array([])
        for i, filter_band in enumerate(self.key):
            dereddened_mag = self.magnitudes[i] - self.av * get_filter_mean_extinction(filter_band)
            expected = np.append(expected, dereddened_mag)
        result = ccm89_deredden_magnitudes(self.key, self.magnitudes, self.av)

        self.assertEqual(expected.tolist(), result.tolist())

    def test_ccm89_dereddens_to_flat_sed(self):
        expected = np.array([0.0, 0.0, 0.0, 0.0, 0.0])
        result = ccm89_deredden_magnitudes(self.key, self.av * self.mean_extinction, self.av)
        self.assertEqual(expected.tolist(), result.tolist())

    def test_build_flux_wl_array_returns_correct_flux_array(self):
        expected = np.array([])
        for i, filter_band in enumerate(self.key):
            flux, effective_wl = mag2flux(filter_band, self.magnitudes[i])
            expected = np.append(expected, flux)
        expected = u.Quantity(expected, flux.unit)

        result_wl_array, result_flux_array = build_flux_wl_array(self.key, 
                                                           self.magnitudes)

        self.assertEqual(expected.tolist(), result_flux_array.tolist())

    def test_build_flux_wl_array_returns_correct_wl_array(self):
        expected = np.array([])
        for i, filter_band in enumerate(self.key):
            flux, effective_wl = mag2flux(filter_band, self.magnitudes[i])
            expected = np.append(expected, effective_wl)
        expected = u.Quantity(expected, effective_wl.unit)

        result_wl_array, result_flux_array = build_flux_wl_array(self.key, 
                                                           self.magnitudes)

        self.assertEqual(expected.tolist(), result_wl_array.tolist())

    def test_integrate_fqbol_returns_correct_flux_value(self):
        expected = np.trapz(self.flux_array, self.eff_wl_array)
        result = integrate_fqbol(self.eff_wl_array, self.flux_array)

        self.assertEqual(expected, result)

    def test_uv_correction_linear(self):
        wavelengths = u.Quantity([0.0, 3660.], u.Angstrom)
        fluxes = u.Quantity([0.0, 5.87565760e-12], u.erg / (u.Angstrom * u.cm**2 * u.s))
        expected = np.trapz(fluxes,wavelengths)
        result = uv_correction_linear(3660. * u.Angstrom,
                                      5.87565760e-12 * u.erg / (u.Angstrom * u.cm**2 * u.s))
        self.assertEqual(expected, result)

class TestBlackbodyIntegration(unittest.TestCase):
        
    def setUp(self):
        self.best_fit_temperature = 4565.75044459
        self.best_fit_angular_radius = 4.59019554625e-09
        self.longest_wavelength = 7980.
        self.shortest_wavelength = 3660.

    def test_ir_correction(self):
        expected = integrate.quad(bb_flux_nounits, self.longest_wavelength,
                                  np.inf,
                                  args=(self.best_fit_temperature,
                                  self.best_fit_angular_radius))
        result = ir_correction(self.best_fit_temperature,
                               self.best_fit_angular_radius,
                               self.longest_wavelength)
        self.assertAlmostEqual(expected[0], result[0])
    
    def test_uv_correction_blackbody(self):
        expected = integrate.quad(bb_flux_nounits, 0, self.shortest_wavelength,
                                  args=(self.best_fit_temperature,
                                        self.best_fit_angular_radius))
        result = uv_correction_blackbody(self.best_fit_temperature,
                                         self.best_fit_angular_radius,
                                         self.shortest_wavelength)
        self.assertAlmostEqual(expected[0], result[0])
