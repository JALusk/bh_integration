import scipy.constants

# Physical constants in cgs units using CODATA2010 measurements
h = (scipy.constants.physical_constants["Planck constant"][0] /
     scipy.constants.erg) # in erg * s
k = (scipy.constants.physical_constants["Boltzmann constant"][0] /
     scipy.constants.erg) # in erg / K
c = (scipy.constants.physical_constants["speed of light in vacuum"][0] *
     1.0E10) # in A / s
