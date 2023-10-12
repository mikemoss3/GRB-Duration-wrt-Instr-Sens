"""
Author: Mike Moss
Contact: mikejmoss3@gmail.com	

Defines the Bayesian block method to calculate the duration of a GRB from a supplied light curve

"""

import numpy as np
from astropy.stats import bayesian_blocks

def bayesian_t_blocks(grb,dur_per=90,ncp_prior=20):
	"""
	Method to extract the duration and photon fluence of a GRB from a supplied light curve using a Bayesian block method. 

	Attributes:
	grb = 		(GRB) 		A grb object 
	dur_per =	(float) 	Percentage of the fluence to calculate the duration for (i.e., T90 corresponds to dur_per = 90)
	"""

	# Bin the light curve 
	bin_edges = bayesian_blocks(t=grb.light_curve['TIME'],x=grb.light_curve['RATE'],sigma=grb.light_curve['UNC'],fitness="measures",ncp_prior=ncp_prior) # Find the T90 and the fluence 

	# Check if any GTI (good time intervals) were found
	if (bin_edges[0] == grb.light_curve['TIME'][0]) and (bin_edges[1] == grb.light_curve['TIME'][-1]):
		# If true, then no GTI's were found
		return 0, 0, 0.
	else:
		# Calculate total duration and start time 
		t_dur_tot = bin_edges[-2] - bin_edges[1]
		t_start_tot = bin_edges[1]

		## Find TXX
		emission_interval = grb.light_curve[np.argmax(t_start_tot<=grb.light_curve['TIME']):np.argmax((t_start_tot+t_dur_tot)<=grb.light_curve['TIME'])]
		# Find the total fluence 
		tot_fluence = np.sum(emission_interval['RATE'])
		# Find the normalized cumulative sum between the total duration 
		cum_sum_fluence = np.cumsum(emission_interval['RATE'])/tot_fluence
		# Find the time interval that encompasses dur_per of the burst fluence
		per_start = ((100 - dur_per)/2)/100
		per_end = 1 - per_start
		t_start =  emission_interval['TIME'][np.argmax(per_start <= cum_sum_fluence)]
		t_end = emission_interval['TIME'][np.argmax(per_end <= cum_sum_fluence)]

		duration = t_end - t_start
		phot_fluence = np.sum(grb.light_curve['RATE'][np.argmax(t_start<=grb.light_curve['TIME']):np.argmax(t_end<=grb.light_curve['TIME'])])

	return duration, t_start, phot_fluence