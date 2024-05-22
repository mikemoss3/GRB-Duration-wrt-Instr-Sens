import numpy as np

def calc_fluence(light_curve, duration, tstart):
	"""
	Method to calculate photon fluence of a light curve

	Attributes:
	---------
	light_curve : nd.array(dtype=[("TIME", float), ("RATE", float), ("UNC", float)])
		Array that stores the light curve 
	duration : float
		Length of the event 
	tstart : float
		What time the event begins
	"""

	dt = light_curve['TIME'][1] - light_curve['TIME'][0]

	# Calculate photon fluence of the light curve within the specified time interval
	fluence = np.sum(light_curve['RATE'][(light_curve['TIME'] > tstart ) & (light_curve['TIME'] < tstart+duration)]) * dt

	# Calculate the 1 second peak flux 
	flux_peak_1s = np.max(light_curve['RATE']) * dt

	return fluence, flux_peak_1s