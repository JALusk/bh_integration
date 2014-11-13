## Description

This code attempts to implement the photometric integration scheme of [Bersten & Hamuy (2009)](http://iopscience.iop.org/0004-637X/701/1/200).

### Method

Bersten and Hamuy calculate bolometric luminosities of three well-observed supernovae in their paper.
In order to do this, they start with photometric observations in multiple filters.
They use the Cardelli et. al. reddening law of $$R_v = 3.1$$ combined with visual absorptions taken from the literature to correct the photometry in each band for reddening.
The magnitudes in each of those filters are then converted to observed fluxes using transmission functions and zero points given in Hamuy (2001).

Once they have fluxes at the effective wavelengths of each photometric band, they use trapezium integration to calculate what they call the quasi-bolometric flux, $$F_{qbol}$$.
When a photometric observation in a certain band isn't available, they interpolate the flux based on the previous and subsequent observations in that band.

In order to go from this quasi-bolometric flux to a true bolometric flux, they need estimates for the flux from the supernova that lies outside the range for which there is photometric coverage.
To do that, they fit a blackbody curve to the observed fluxes.
They do, however, mention that they exclude from these fits observations in the U-band which fall below the blackbody model.
But how do they know a point falls below the model without first including the point in the fit?
My guess is that they start by fitting the reddest fluxes, and including bluer bands one at a time, with some threshold for "if you fall below the BB curve that fits all the data redward of you, then you are excluded from the fit."

Corrections to the quasi-bolometric flux:

To correct for missing flux in the IR, the flux redward of the reddest band observed is approximated using the integral of the best-fitting blackbody function from the effective wavelength of the reddest observation to infinity.
For their well-observed supernovae, this correction $$F_{IR}$$  increases with time as the SED shifts to the red, but always remains below 7%

To correct for missing flux in the UV is more complicated.
There is strong line blanketing in the blue bands of a supernova which will cause the SED to depart from a black body.
This is why the U and B band observations may fall below a BB curve and be excluded from the fit as described above.

Bersten and Hamuy use a number of different techniques to correct for the flux on the UV side of the SED.
If the U-band flux falls on the best-fitting blackbody curve, then the UV correction is calculated by integrating the blackbody from the effective wavelength of the U-band to zero.
If the U-band flux falls below the best-fitting blackbody durve, then the UV correction is calculated by integrating the blackbody from the effective wavelength of the U-band along a straight line to zero flux at $$\lambda = 2000\r{A}$$.
This is done based on the behavior of the SED of supernova atmosphere models.

In order to determine the reliability of these UV corrections, the authors compare them to the atmosphere models of Eastman et. al. (1996) and Dessart & Hillier (2005).
Those models show that the UV correction at very early times is higher than that calculated from the method described above.
They use the $$F_{UV}$$ from the atmosphere models when the $$B-V$$ color of the supernova is $$< -0.04$$.

### Current Status

Tests passing.

All the major pieces of the code (save for those noted in the ToDo section) are written and are passing code tests.

There are two functions which return the blackbody flux. This is because the integration functions that I am using don't play nice with astropy quantities. To save the headache of maintaining two copes of the same function, bb_flux() returns the astropy quantity, and bb_flux_nounits() is a dummy wrapper function that just calls bb_flux() and then strips off the units, returning a float. 
I chose to handle the problem in this way so that there is only one function that truly calculates the blackbody flux, thus making testing easier and lessening the possibility that bugs get introduced that affect one flux calculation and not the other.
This is the function called by the scipy.integrate.quad() routine during the calculation of the IR and UV flux corrections.

### ToDo:

    1) De-reddening routine
    2) Logic that picks which UV correction to use based on whether or not U lies below the BB fit
    3) Docstrings for all functions

