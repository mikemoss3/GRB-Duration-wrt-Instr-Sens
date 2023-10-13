"""
Author: Mike Moss
Contact: mikejmoss3@gmail.com	

Defines all the functions necessary to simulate an observation of a GRB using an input template, designated instrument response, and observing conditions 

"""

import numpy as np
from packages.class_GRB import GRB
from util_packages.package_det_ang_dependence import find_pcode, find_inc_ang


def simulate_observation(template_grb, imx, imy, ndets, resp_mat, z_p=0, sim_triggers=False,ndet_max=32768,bgd_rate_per_det=0.3):
	"""
	Method to complete a simulation of a synthetic observation based on the input source frame GRB template and the desired observing conditions

	Attributes:
	template_grb : GRB 
		GRB class object that holds the source frame information of the template GRB
	imx, imy : 			float, float 
		The x and y position of the GRB on the detector plane
	resp_mat : RSP
		Response matrix to convolve the template spectrum with
	ndets : int
		Number of detectors enabled during the synthetic observation 
	z : float 
		Redshift of synthetic GRB
	sim_triggers : boolean
		Whether or not to simulate the Swift/BAT trigger algorithms or not
	ndet_max : int
		Maximum number of detectors on the detector plane (for Swift/BAT ndet_max = 32,768)
	bdg_rate_per_det : float 
		Background level to be added to the synthetic light curve
	"""

	# Initialize synth_GRB
	synth_GRB = template_grb.copy()
	synth_GRB.imx, synth_GRB.imy = imx, imy
	synth_GRB.z = z_p

	# Apply distance corrections to template GRB light curve to create synthetic GRB light cure.
	# This assumes that the template light curve is in the source frame
	synth_GRB.move_to_new_frame(z_o=template_grb.z, z_p=z_p)

	# Apply observing condition corrections (e.g., NDETS)
	det_frac = ndets / ndet_max # Current number of enabled detectors divided by the maximum number of possible detectors

	# Fold GRB through instrument response (RSP selected based on position on the detector plane)
	folded_spec = resp_mat.fold_spec(template_grb.specfunc)
	rate_15_350keV = band_rate(folded_spec,15,350) * det_frac
	synth_GRB.light_curve['RATE'] = synth_GRB.light_curve['RATE']*rate_15_350keV
	synth_GRB.light_curve['UNC'] = synth_GRB.light_curve['UNC']*rate_15_350keV

	# If we are testing the trigger algorithm:
		# Modulate the light curve by the folded spectrum normalization for each energy band 
		# Calculate the fraction of the quadrant exposure 

	# Add background to light curve 
	bgd_rate = bgd_rate_per_det * ndets 
	synth_GRB.light_curve['RATE'] += bgd_rate

	# Apply fluctuations 
	# synth_GRB.light_curve['RATE'] = np.random.normal(loc=synth_GRB.light_curve['RATE'],scale=synth_GRB.light_curve['UNC'])
	synth_GRB.light_curve['RATE'] = np.random.normal(loc=synth_GRB.light_curve['RATE'],scale=np.sqrt(synth_GRB.light_curve['RATE']))

	# Apply mask-weighting to light curve (both the rate and uncertainty)
	synth_GRB.light_curve = apply_mask_weighting(synth_GRB.light_curve,imx,imy,ndets,bgd_rate)

	return synth_GRB

def band_rate(spectrum,emin,emax):
	"""
	Method to calculate the rate by taking the sum of the spectrum across a specified energy band
	"""

	return np.sum(spectrum['RATE'][np.argmax(spectrum['ENERGY']>=emin):np.argmax(spectrum['ENERGY']>=emax)])


def apply_mask_weighting(light_curve,imx,imy,ndets,bgd_rate):
	"""
	Method to apply mask weighting to a light curve assuming a flat background.
	Mask-weighted means:
		1. Background subtraction
		2. Per detector
		3. Per illuminated detector (partial coding fraction)
		4. Fraction of detector illuminated (mask correction)
		5. On axis equivalent (effective area correction for off-axis bursts)
	"""

	# From imx and imy, find pcode and the angle of incidence
	pcode = find_pcode(imx,imy)
	angle_inc = find_inc_ang(imx,imy) # rad

	# Total mask-weighting correction
	correction = np.cos(angle_inc*np.pi/180)*pcode*ndets*fraction_correction(pcode) # total correction factor

	# Rough calculation of the background standard deviation 
	# stdev_backgroud = np.std(light_curve["RATE"][0:20])
	stdev_backgroud = np.sqrt(np.mean(light_curve['RATE'][0:20]))

	# Calculate the mask-weighted RATE column
	light_curve['RATE'] = (light_curve["RATE"] - bgd_rate)/correction
	# Use error propagation to calculate the uncertainty in the RATE for the mask-weight light curve
	light_curve['UNC'] = np.sqrt( np.power(light_curve['UNC']/correction,2.)+np.power(stdev_backgroud/correction,2.))

	return light_curve

def fraction_correction(pcode):
	"""
	Method that calculates and returns a correction fraction that was found to be needed for off-axis bursts. 
	This factor is needed to correct for the FFT convolution that is used for Swift/BAT
	This correction was empirically fit with a quadratic function, which is how the parameter values in this method were determined.
	"""
	a=1.1205830634986802
	b=-1.2137924102819533
	c=0.6178450561688628
	return a + (b*pcode) + (c*np.power(pcode,2))
