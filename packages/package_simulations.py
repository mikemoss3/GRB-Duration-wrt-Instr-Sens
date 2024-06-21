"""
Author: Mike Moss
Contact: mikejmoss3@gmail.com	

Defines all the functions necessary to simulate an observation of a GRB using an input template, designated instrument response, and observing conditions 

"""

import numpy as np
from scipy.stats import rv_discrete
from packages.class_GRB import GRB
from packages.class_RSP import ResponseMatrix
from util_packages.package_det_ang_dependence import find_pcode, find_inc_ang, fraction_correction


def simulate_observation(template_grb, synth_grb, imx, imy, ndets, resp_mat,
	z_p=0, sim_triggers=False,
	ndet_max=32768, bgd_rate_per_det=0.3, band_rate_min=14, band_rate_max=350):
	"""
	Method to complete a simulation of a synthetic observation based on the input source frame GRB template and the desired observing conditions

	Attributes:
	------------------------
	template_grb : GRB 
		GRB class object that holds the source frame information of the template GRB
	synth_grb : GRB 
		GRB class object that will hold the simulated light curve
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
	"""

	import matplotlib.pyplot as plt
	from astropy.io import fits


	# Initialize synth_GRB
	synth_grb.imx, synth_grb.imy = imx, imy
	synth_grb.z = z_p

	# Apply distance corrections to GRB light curve and spectrum
	synth_grb.move_to_new_frame(z_o=template_grb.z, z_p=z_p)

	# Calculate the fraction of the detectors currently enabled 
	det_frac = ndets / ndet_max # Current number of enabled detectors divided by the maximum number of possible detectors

	folded_spec = resp_mat.fold_spec(synth_grb.specfunc)  # counts / sec / keV
	rate_in_band = band_rate(folded_spec, band_rate_min, band_rate_max) * det_frac # counts / sec

	# Using the total count rate from the spectrum and the relative flux level of the light curve, make a new light curve
	# The synthetic GRB light curve technically has units of counts / sec / cm^2, but we are only using it as a template for relative flux values. 
	# The actual flux is set by the band rate, which is in units of counts / sec 
	synth_grb.light_curve['RATE'] = synth_grb.light_curve['RATE'] * rate_in_band # counts / sec

	# If we are testing the trigger algorithm:
		# Modulate the light curve by the folded spectrum normalization for each energy band 
		# Calculate the fraction of the quadrant exposure 

	# Apply mask-weighting to light curve (both the rate and uncertainty)
	synth_grb.light_curve = apply_mask_weighting(synth_grb.light_curve, imx, imy, ndets, 0) # background-subtracted counts / sec / det

	# The length is determined by the light curve time bin size
	sim_lc_length = int( (2*template_grb.t_buffer/template_grb.dt) + len(synth_grb.light_curve) )

	# Initialize an empty light curve to hold the simulated light curve 
	empty = np.zeros(shape=sim_lc_length, dtype=[('TIME', float), ('RATE',float), ('UNC',float)])
	# Fill the time axis from synth_grb-buffer to synth_grb+buffer with correct time bin sizes 
	empty['TIME'] = np.arange(
		start=synth_grb.light_curve['TIME'][0]-template_grb.t_buffer, 
		stop= synth_grb.light_curve['TIME'][-1]+template_grb.t_buffer+template_grb.dt, 
		step=template_grb.dt
		)[:len(empty)]

	# Pull a random background variance from the distribution created observed values
	variance = rand_background_variance()
	# Correct for time-bin size
	variance /= np.sqrt(template_grb.dt)

	empty['RATE'] = np.random.normal( loc=np.zeros(shape=len(empty)), scale=variance)
	# Set the uncertainty of the count rate to the standard deviation. 
	empty['UNC'] = np.ones(shape=len(empty))*variance

	len_sim = len(synth_grb.light_curve['RATE'])
	argstart = np.argmax(empty['TIME']>=synth_grb.light_curve['TIME'][0])
	empty[argstart: argstart+len_sim]['RATE'] += synth_grb.light_curve['RATE']

	synth_grb.light_curve = empty

	return synth_grb

def band_rate(spectrum,emin,emax):
	"""
	Method to calculate the rate by taking the sum of the spectrum across a specified energy band
	"""

	return np.sum(spectrum['RATE'][np.argmax(spectrum['ENERGY']>=emin):np.argmax(spectrum['ENERGY']>=emax)])


def apply_mask_weighting(light_curve, imx, imy, ndets, bgd_rate):
	"""
	Method to apply mask weighting to a light curve assuming a flat background.
	Mask-weighted means:
		1. Background subtraction
		2. Per detector
		3. Per illuminated detector (partial coding fraction)
		4. Fraction of detector illuminated (mask correction)
		5. On axis equivalent (effective area correction for off-axis bursts)
	"""

	# Rough calculation of the background standard deviation 
	stdev_backgroud = np.sqrt(bgd_rate)

	# From imx and imy, find pcode and the angle of incidence
	pcode = find_pcode(imx, imy)
	angle_inc = find_inc_ang(imx,imy) # rad

	# Total mask-weighting correction
	correction = np.cos(angle_inc)*pcode*ndets*fraction_correction(imx, imy) # total correction factor

	if pcode == 0:
		# Source was not in the field of view
		light_curve['UNC'] *= 0
		light_curve['RATE'] *=0 # counts / sec / dets
		return light_curve

	# Use error propagation to calculate the uncertainty in the RATE for the mask-weight light curve
	light_curve['UNC'] = np.sqrt( np.power(np.sqrt(np.abs(light_curve['RATE']))/correction,2.)+np.power(stdev_backgroud/correction,2.))
	# Calculate the mask-weighted RATE column
	light_curve['RATE'] = (light_curve["RATE"] - bgd_rate)/correction # background-subtracted counts / sec / dets


	return light_curve

def add_background(light_curve, bgd_rate, buffer, dt):
	"""
	Method to add a buffer interval and a flat background to a given light curve

	Attributes:
	------------------------
	light_curve : np.ndarray([("TIME",float), ("RATE",float), ("UNC",float)])
		Light curve array
	bgd_rate : float
		Background rate (counts / sec) to be added to light curve 
	buffer : float
		Duration (sec) of the background interval to be added to either side of the existing light curve
	dt : float
		time bin size
	"""

	old_lc_size = len(light_curve)
	old_time_start = light_curve['TIME'][0]
	old_time_end = light_curve['TIME'][-1]
	buffer_size = int(np.floor(buffer/dt))

	# Add buffer interval to light curve 
	buffer_interval = np.zeros(shape=buffer_size, dtype=[('TIME',float), ('RATE',float), ('UNC',float)])
	light_curve = np.concatenate( (buffer_interval, light_curve, buffer_interval), axis=0) 

	# Fill time axis value
	light_curve['TIME'][:buffer_size] = np.arange(old_time_start - buffer, old_time_start, step=dt)[:buffer_size]
	light_curve['TIME'][buffer_size+old_lc_size:] = np.arange(old_time_end+dt, old_time_end+buffer+dt, step=dt)[:buffer_size]

	# Add flat background
	light_curve['RATE'] += bgd_rate

	return light_curve

def fit_function(t,fm,tm,r,d):
	"""
	FRED shaped light curve based on Equation 22 from Kocevski et al. 2003, basd on power law rise and exponential decay

	Attributes:
	t = time since trigger
	fm = flux at peak of the pulse (fm = F(tm))
	tm = t_max or the peak time of the pulse 
	r = rise constant
	d = decay constant 

	"""

	flux = fm*np.power(t/tm,r)*np.power( (d/(d+r)) + ((r/(d+r))*np.power(t/tm,r+1)) ,-(r+d)/(r+1))

	return flux

def rand_background_variance():
	"""
	Method that return a randomly selected from a distribution created from the measured background variances observed by Swift/BAT.

	The distribution function is created from the FRED function described in Kocevski 2003. 
	The parameter values were found in a separate fit.  
	"""

	# There are only anomolous variances found outside of these cuts
	cut_min = 0.02
	cut_max = 0.25

	x_range = np.linspace(cut_min, cut_max)

	# parameters = 0.0869, 0.0686, 9.929, 3.066 
	parameters = [0.08694067, 0.06859141, 9.92930295, 3.06620674]
	"""
	Parameters:
	fm = flux at peak of the pulse (fm = F(tm))
	tm = t_max or the peak time of the pulse 
	r = rise constant
	d = decay constant
	"""

	# Create distribution
	distrib = rv_discrete(a=cut_min, b=cut_max, values=(x_range, fit_function(x_range, *parameters)) ) 
	# Select random value
	variance = distrib.rvs(size=1)

	return variance