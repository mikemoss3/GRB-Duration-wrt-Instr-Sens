"""
Author: Mike Moss
Contact: mikejmoss3@gmail.com	

Defines the Bayesian block method to calculate the duration of a GRB from a supplied light curve

"""

import numpy as np
import scipy 
from astropy.stats import bayesian_blocks
# Import the time library
import time

def bayesian_t_blocks(light_curve, dur_per=90, ncp_prior=6):
	"""
	Method to extract the duration and photon fluence of a GRB from a supplied light curve using a Bayesian block method. 

	Attributes:
	grb = 		(GRB) 		A grb object 
	dur_per =	(float) 	Percentage of the fluence to calculate the duration for (i.e., T90 corresponds to dur_per = 90)
	"""

	# Astropy bayesian block algorithm -- it is much safer to use since it handles exceptions better, however its an order of magnitude
	# bin_edges = bayesian_blocks(t=light_curve['TIME'], x=light_curve['RATE'], sigma=light_curve['UNC'], fitness="measures", ncp_prior=ncp_prior) # Find the T90 and the fluence 
	# condensed bayesian block algorithm taken from astropy (faster but more unsafe)
	bin_edges = custom_bb(light_curve, ncp_prior)

	# Check if any GTI (good time intervals) were found
	if len(bin_edges) <= 3:
		# If true, then no GTI's were found		
		return 0., 0., 0.
	else:
		# Calculate total duration and start time 
		t_dur_tot = bin_edges[-2] - bin_edges[1]
		t_start_tot = bin_edges[1]

		## Find TXX
		emission_interval = light_curve[np.argmax(t_start_tot<=light_curve['TIME']):np.argmax((t_start_tot+t_dur_tot)<=light_curve['TIME'])]
		if len(emission_interval) == 0:
			# Then no Bayesian blocks were found.
			return 0, 0, 0
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
		phot_fluence = np.sum(light_curve['RATE'][np.argmax(t_start<=light_curve['TIME']):np.argmax(t_end<=light_curve['TIME'])])

	return duration, t_start, phot_fluence

def custom_bb(light_curve, ncp_prior):

	t=light_curve['TIME']
	x=light_curve['RATE']
	sigma=light_curve['UNC']

	ak_raw = np.ones_like(x) / sigma**2
	bk_raw = x / sigma**2

	edges = np.concatenate([t[:1], 0.5 * (t[1:] + t[:-1]), t[-1:]])

	# arrays to store the best configuration
	N = len(t)
	best = np.zeros(N, dtype=float)
	last = np.zeros(N, dtype=int)

	# ----------------------------------------------------------------
	# Start with first data cell; add one cell at each iteration
	# ----------------------------------------------------------------
	for R in range(N):
		# a_k: eq. 31
		a_k = 0.5 * np.cumsum(ak_raw[: (R + 1)][::-1])[::-1]

		# b_k: eq. 32
		b_k = -np.cumsum(bk_raw[: (R + 1)][::-1])[::-1]

		# evaluate fitness function
		fit_vec = fitness(a_k, b_k)

		A_R = fit_vec - ncp_prior
		A_R[1:] += best[:R]

		i_max = np.argmax(A_R)
		last[R] = i_max
		best[R] = A_R[i_max]

	# ----------------------------------------------------------------
	# Now find changepoints by iteratively peeling off the last block
	# ----------------------------------------------------------------
	change_points = np.zeros(N, dtype=int)
	i_cp = N
	ind = N
	while i_cp > 0:
	    i_cp -= 1
	    change_points[i_cp] = ind
	    if ind == 0:
	        break
	    ind = last[ind - 1]
	if i_cp == 0:
	    change_points[i_cp] = 0
	change_points = change_points[i_cp:]

	return edges[change_points]


def fitness(a_k, b_k):
	# eq. 41 from Scargle 2013
	return (b_k * b_k) / (4 * a_k)