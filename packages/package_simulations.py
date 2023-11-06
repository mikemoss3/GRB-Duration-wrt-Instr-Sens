"""
Author: Mike Moss
Contact: mikejmoss3@gmail.com	

Defines all the functions necessary to simulate an observation of a GRB using an input template, designated instrument response, and observing conditions 

"""

import numpy as np
from packages.class_GRB import GRB
from packages.class_RSP import ResponseMatrix
from util_packages.package_det_ang_dependence import find_pcode, find_inc_ang, fraction_correction


def simulate_observation(template_grb, imx, imy, ndets, 
	resp_mat=None, z_p=0, sim_triggers=False,
	ndet_max=32768, bgd_rate_per_det=0.3, area_per_det = 0.16, band_rate_min=15,band_rate_max=350):
	"""
	Method to complete a simulation of a synthetic observation based on the input source frame GRB template and the desired observing conditions

	Attributes:
	template_grb : GRB 
		GRB class object that holds the source frame information of the template GRB
	imx, imy : 	float, float 
		The x and y position of the GRB on the detector plane
	resp_mat : RSP
		Response matrix to convolve the template spectrum with. If no response matrix is given, a Swift/BAT response matrix is assumed from the given imx, imy
	ndets : int
		Number of detectors enabled during the synthetic observation 
	z_p : float 
		Redshift of synthetic GRB
	sim_triggers : boolean
		Whether or not to simulate the Swift/BAT trigger algorithms or not
	ndet_max : int
		Maximum number of detectors on the detector plane (for Swift/BAT ndet_max = 32,768)
	bdg_rate_per_det : float 
		Background level to be added to the synthetic light curve
	area_per_det : float
		Area of a single detector element in cm^2
	"""

	# Initialize synth_GRB
	synth_GRB = template_grb.copy() # Copies light curve and spectrum
	synth_GRB.imx, synth_GRB.imy = imx, imy
	synth_GRB.z = z_p


	# Apply distance corrections to GRB light curve and spectrum
	synth_GRB.move_to_new_frame(z_o=template_grb.z, z_p=z_p)

	# Calculate the fraction of the detectors currently enabled 
	det_frac = ndets / ndet_max # Current number of enabled detectors divided by the maximum number of possible detectors

	# Fold GRB through instrument response (RSP selected based on position on the detector plane)
	if resp_mat is None:
		resp_mat = ResponseMatrix()
		resp_mat.load_SwiftBAT_resp(imx,imy)


	folded_spec = resp_mat.fold_spec(synth_GRB.specfunc)  # counts / sec / keV
	rate_in_band = band_rate(folded_spec, band_rate_min, band_rate_max) * det_frac # counts / sec 

	# Using the total count rate from the spectrum and the relative flux level of the light curve, make a new light curve
	# The synthetic GRB light curve technically has units of counts / sec / cm^2, but we are only using it as a template for relative flux values. 
	# The actual flux is set by the band rate, which is in units of counts / sec 
	synth_GRB.light_curve['RATE'] = synth_GRB.light_curve['RATE'] * rate_in_band # counts / sec
	synth_GRB.light_curve['UNC'] = synth_GRB.light_curve['UNC'] * rate_in_band

	# If we are testing the trigger algorithm:
		# Modulate the light curve by the folded spectrum normalization for each energy band 
		# Calculate the fraction of the quadrant exposure 

	# Add background to light curve 
	bgd_rate = bgd_rate_per_det * ndets # counts / sec
	synth_GRB.light_curve['RATE'] += bgd_rate # counts / sec

	# Apply fluctuations 
	synth_GRB.light_curve['RATE'] = np.random.normal(loc=synth_GRB.light_curve['RATE'],scale=np.sqrt(synth_GRB.light_curve['RATE'])) # counts / sec

	# Apply mask-weighting to light curve (both the rate and uncertainty)
	synth_GRB.light_curve = apply_mask_weighting(synth_GRB.light_curve,imx,imy,ndets,bgd_rate) # background-subtracted counts / sec / det

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
	correction = np.cos(angle_inc)*pcode*ndets*fraction_correction(imx, imy) # total correction factor

	# Rough calculation of the background standard deviation 
	stdev_backgroud = np.sqrt(np.mean(light_curve['RATE'][0:20]))

	# Use error propagation to calculate the uncertainty in the RATE for the mask-weight light curve
	light_curve['UNC'] = np.sqrt( (light_curve['RATE']/np.power(correction,2.))+np.power(stdev_backgroud/correction,2.))
	# Calculate the mask-weighted RATE column
	light_curve['RATE'] = (light_curve["RATE"] - bgd_rate)/correction # background-subtracted counts / sec / dets

	return light_curve
