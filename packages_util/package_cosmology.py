"""
Author: Mike Moss
Contact: mikejmoss3@gmail.com	
Last edited: 2023-10-12

This file defines functions that deal with cosmological corrections.

"""

import numpy as np
import scipy.integrate as integrate 
import packages_util.globalconstants as gc


def lum_dis(z: float):
	""" 
	Caclulate luminosity distance for redshift z
	"""
	lum_dis_Mpc = ((1+z)*gc.c/(gc.H0) ) * integrate.quad(lambda zi: 1/np.sqrt( ((gc.omega_m*np.power(1+zi,3) )+gc.omega_lam) ),0,z)[0]
	lum_dis_cm = lum_dis_Mpc * 3.086e24 # Mpc -> cm
	return lum_dis_cm

def k_corr(func, params, z, en_band_min, en_band_max):
	""" 
	Calculates the bolumetric k-correction using a specified function form at a particular redshft
	
	Attributes:
	func = spectral function
	params = function parameters
	z = redshift

	""" 

	numerator = integrate.quad(lambda en: en*func(en, params), gc.bol_lum[0],gc.bol_lum[1])[0]
	denominator = integrate.quad(lambda en: en*func(en, params), en_band_min*(1+z), en_band_max*(1+z))[0]
	
	return numerator/denominator