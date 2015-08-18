import numpy as np
from astropy import constants as const
from astropy import units as u
# from planck import *
import my_funcs
from scipy import integrate
h = const.h.to('erg s').value
c = const.c.to('cm/s').value
k_b = const.k_B.to('erg/(K)').value
ryd = const.Ryd.to('1/(cm)').value # h*c*ryd = 13.6ev
sigma_sb = const.sigma_sb.to("erg/(s K4 cm2)").value

#-- write planck function in dimensionless form x = h\nu/kT
#-- Also write it so that it doesn't blow up at x --> \infinity
def Bnu(x,T):
  return 2.0/(h**3*c**2)*(k_b*T)**4*x**3*np.exp(-x)/(1.0-np.exp(-x))

def test_integral(temp):
  # def planck_function_nounits(wavelength, temperature):
  #   ans = planck_function(wavelength, temperature)
  #   return ans.value
  T = temp 
  nenner,err = integrate.quad(Bnu,0.,np.inf,args=(T,))   # this is just sigma/pi T**4 
  result1,err = integrate.quad(my_funcs.planck_wl,100.,np.inf,args=(T,))
  result2,err = integrate.quad(my_funcs.planck_nu,1.e13,1.e17,args=(T,))
  expected = const.sigma_sb.cgs * (T*u.K)**4/np.pi
  return nenner,result1,result2,expected.value
