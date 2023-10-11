"""
Author: Mike Moss
Contact: mikejmoss3@gmail.com	
Last edited: 2023-10-12

This file defines functions that deal with cosmological corrections.

"""

import numpy as np
import scipy.integrate as integrate 
import util_packages.globalconstants as gc


def lum_dis(z: float):
	""" 
	Caclulate luminosity distance for redshift z
	"""
	lum_dis_Mpc = ((1+z)*gc.c/(gc.H0) ) * integrate.quad(lambda zi: 1/np.sqrt( ((gc.omega_m*np.power(1+zi,3) )+gc.omega_lam) ),0,z)[0]
	lum_dis_cm = lum_dis_Mpc * 3.086e24 # Mpc -> cm
	return lum_dis_cm

def k_corr(specfunc, z, emin, emax):
	""" 
	Calculates the bolumetric k-correction using a specified function form at a particular redshft. See Bloom, Frail, and Sari 2001.
	
	Attributes:
	func = spectral function
	params = function parameters
	z = redshift

	""" 

	# Create restframe spectrum by correcting E_peak or temperature by the redshift (if spectral function has a peak energy or temperature)
	rs_specfunc = specfunc.copy()
	for i, (key, val) in enumerate(rs_specfunc.params.items()):
		if key == "ep":
			rs_specfunc.params[key] *= (1+z)
		if key == "temp":
			rs_specfunc.params[key] *= (1+z)

	# Evaluate bolometric spectrum in the rest frame of the source 
	numerator = integrate.quad(lambda en: en*rs_specfunc(en), gc.bol_lum[0]/(1+z),gc.bol_lum[1]/(1+z))[0]
	# Evaluate spectrum within the defined band pass in the observer frame
	denominator = integrate.quad(lambda en: en*specfunc(en), emin, emax)[0]
	
	return numerator/denominator