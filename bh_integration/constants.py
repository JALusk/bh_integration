import scipy.constants

# Physical constants in cgs units using CODATA2010 measurements
h = (scipy.constants.physical_constants["Planck constant"][0] /
     scipy.constants.erg)
k = (scipy.constants.physical_constants["Boltzmann constant"][0] /
     scipy.constants.erg)
c = (scipy.constants.physical_constants["speed of light in vacuum"][0] *
     100.0)
