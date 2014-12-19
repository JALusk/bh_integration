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
        self.magnitudes = np.array([7.129, 5.554, 4.383, 3.917, 3.794])
        self.uncertainties = np.array([0.02, 0.02, 0.02, 0.02, 0.02])
        self.av = 0.435
        self.flux_array = u.Quantity([5.87565760e-12, 3.79417256e-11, 
                                      6.40953095e-11, 5.90280490e-11, 
                                      3.41930932e-11], 
                                      u.erg / (u.Angstrom * u.cm**2 * u.s))
        self.flux_uncertainties = u.Quantity([1.0e-12, 1.0e-12, 
                                      1.0e-12, 1.0e-12, 
                                      1.0e-12], 
                                      u.erg / (u.Angstrom * u.cm**2 * u.s))

        self.eff_wl_array = u.Quantity([3660., 4380., 5450., 6410., 7980.],
                                        u.Angstrom)
    
    def test_build_flux_wl_array_returns_correct_flux_array(self):
        expected = np.array([])
        for i, filter_band in enumerate(self.key):
            flux, uncertainty, effective_wl = mag2flux(filter_band, self.magnitudes[i], 
                                                       self.uncertainties[i])
            expected = np.append(expected, flux)
        expected = u.Quantity(expected, flux.unit)

        result_wl_array, result_flux_array, result_flux_undertainty_array = build_flux_wl_array(self.key, 
                                                           self.magnitudes,
                                                           self.uncertainties)

        self.assertEqual(expected.value.tolist(), result_flux_array.value.tolist())

    def test_build_flux_wl_array_returns_correct_flux_uncertainty_array(self):
        expected = np.array([])
        for i, filter_band in enumerate(self.key):
            flux, uncertainty, effective_wl = mag2flux(filter_band, self.magnitudes[i], 
                                                       self.uncertainties[i])
            expected = np.append(expected, uncertainty)
        expected = u.Quantity(expected, flux.unit)

        result_wl_array, result_flux_array, result_flux_uncertainty_array = build_flux_wl_array(self.key, 
                                                           self.magnitudes,
                                                           self.uncertainties)

        self.assertEqual(expected.value.tolist(), result_flux_uncertainty_array.value.tolist())

    def test_build_flux_wl_array_returns_correct_wl_array(self):
        expected = np.array([])
        for i, filter_band in enumerate(self.key):
            flux, uncertainty, effective_wl = mag2flux(filter_band, self.magnitudes[i],
                                                       self.uncertainties[i])
            expected = np.append(expected, effective_wl)
        expected = u.Quantity(expected, effective_wl.unit)

        result_wl_array, result_flux_array, result_flux_uncertainty_array = build_flux_wl_array(self.key, 
                                                           self.magnitudes,
                                                           self.uncertainties)

        self.assertEqual(expected.value.tolist(), result_wl_array.value.tolist())

    def test_integrate_fqbol_returns_correct_flux_value(self):
        expected = np.trapz(self.flux_array, self.eff_wl_array)
        result, uncertainty = integrate_fqbol(self.eff_wl_array, self.flux_array, self.flux_uncertainties)

        self.assertEqual(expected, result)

    def test_integrate_fqbol_returns_correct_flux_uncertainty(self):
        quad_terms = np.array([])

        for i, uncertainty in enumerate(self.flux_uncertainties):
            if i == 0:
                term = 0.5 * (self.eff_wl_array[i+1] - self.eff_wl_array[i]) * uncertainty
                quad_terms = np.append(quad_terms, term)
            elif i == len(self.flux_uncertainties) - 1:
                term = 0.5 * (self.eff_wl_array[i] - self.eff_wl_array[i-1]) * uncertainty
                quad_terms = np.append(quad_terms, term)
            else:
                term = 0.5 * (self.eff_wl_array[i+1] - self.eff_wl_array[i-1]) * uncertainty
                quad_terms = np.append(quad_terms, term)

        expected = np.sqrt(np.sum(x*x for x in quad_terms))
        result, uncertainty = integrate_fqbol(self.eff_wl_array, self.flux_array, self.flux_uncertainties)

        self.assertEqual(expected, uncertainty.value)
    def test_uv_correction_linear(self):
        wavelengths = u.Quantity([0.0, 3660.], u.Angstrom)
        fluxes = u.Quantity([0.0, 5.87565760e-12], u.erg / (u.Angstrom * u.cm**2 * u.s))
        expected = np.trapz(fluxes,wavelengths)
        result = uv_correction_linear(3660. * u.Angstrom,
                                      5.87565760e-12 * u.erg / (u.Angstrom * u.cm**2 * u.s))
        self.assertEqual(expected, result)

    def test_calculate_fbol(self):
        wavelengths, fluxes, flux_uncertainties = build_flux_wl_array(self.key, self.magnitudes, self.uncertainties)
        fluxes_unred = fluxes * extinction.reddening(wavelengths, self.av, model='ccm89')

        shortest_wl = np.amin(wavelengths)
        longest_wl = np.amax(wavelengths)
        
        fqbol, fqbol_uncertainty = integrate_fqbol(wavelengths, fluxes_unred, flux_uncertainties)

        temperature, angular_radius, chisq = bb_fit_parameters(wavelengths, 
                                                        fluxes_unred, flux_uncertainties)
        
        ir_corr = ir_correction(temperature.value, angular_radius, longest_wl.value)[0]
        uv_corr = uv_correction_blackbody(temperature.value, angular_radius, 
                                          shortest_wl.value)[0]

        expected = fqbol.value + ir_corr + uv_corr

        result, uncertainty = calculate_fbol(self.key, self.magnitudes, self.uncertainties, self.av)

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
    
    def test_ir_correction_uncertainty(self):
        expected = integrate.quad(bb_flux_nounits, self.longest_wavelength,
                                  np.inf,
                                  args=(self.best_fit_temperature,
                                  self.best_fit_angular_radius))
        result = ir_correction(self.best_fit_temperature,
                               self.best_fit_angular_radius,
                               self.longest_wavelength)
        self.assertAlmostEqual(expected[1], result[1])

    def test_uv_correction_blackbody(self):
        expected = integrate.quad(bb_flux_nounits, 0, self.shortest_wavelength,
                                  args=(self.best_fit_temperature,
                                        self.best_fit_angular_radius))
        result = uv_correction_blackbody(self.best_fit_temperature,
                                         self.best_fit_angular_radius,
                                         self.shortest_wavelength)
        self.assertAlmostEqual(expected[0], result[0])

    def test_uv_correction_blackbody_uncertainty(self):
        expected = integrate.quad(bb_flux_nounits, 0, self.shortest_wavelength,
                                  args=(self.best_fit_temperature,
                                        self.best_fit_angular_radius))
        result = uv_correction_blackbody(self.best_fit_temperature,
                                         self.best_fit_angular_radius,
                                         self.shortest_wavelength)
        self.assertAlmostEqual(expected[1], result[1])
