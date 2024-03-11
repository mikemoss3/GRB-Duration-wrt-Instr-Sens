"""
Author: Mike Moss
Contact: mikejmoss3@gmail.com	

Defines the Bayesian block method to calculate the duration of a GRB from a supplied light curve

"""

import numpy as np
import scipy 
from astropy.stats import bayesian_blocks


def bayesian_t_blocks(grb, dur_per=90, ncp_prior=6):
	"""
	Method to extract the duration and photon fluence of a GRB from a supplied light curve using a Bayesian block method. 

	Attributes:
	grb = 		(GRB) 		A grb object 
	dur_per =	(float) 	Percentage of the fluence to calculate the duration for (i.e., T90 corresponds to dur_per = 90)
	"""

	# Bin the light curve 
	bin_edges = bayesian_blocks(t=grb.light_curve['TIME'], x=grb.light_curve['RATE'], sigma=grb.light_curve['UNC'], fitness="measures", ncp_prior=ncp_prior) # Find the T90 and the fluence 

	# Check if any GTI (good time intervals) were found
	if len(bin_edges) <= 3:
		# If true, then no GTI's were found

		# Set duration information for the GRB object 
		grb.set_duration(0., 0., 0., dur_per, ncp_prior)
		
		return 0., 0., 0.
	else:
		# Calculate total duration and start time 
		t_dur_tot = bin_edges[-2] - bin_edges[1]
		t_start_tot = bin_edges[1]

		## Find TXX
		emission_interval = grb.light_curve[np.argmax(t_start_tot<=grb.light_curve['TIME']):np.argmax((t_start_tot+t_dur_tot)<=grb.light_curve['TIME'])]
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
		phot_fluence = np.sum(grb.light_curve['RATE'][np.argmax(t_start<=grb.light_curve['TIME']):np.argmax(t_end<=grb.light_curve['TIME'])])

		# Set duration information for the GRB object 
		grb.set_duration(duration, t_start, phot_fluence, dur_per, ncp_prior)


	return duration, t_start, phot_fluence


def custom_bblocks(light_curve, ncp_prior=6, coalescefrac = 0.05):
	"""
	Bayesian block implementation based on heasoft tool battblocks
	"""

	nlag = 0
	timedel = light_curve['TIME'][1] - light_curve['TIME'][0]
	nrows = len(light_curve)

	# Convert light curve to cells with size in units of timebins. 
	# Since we are using constant time bins, each cell has a size of one
	cell_curve = np.ones(shape=len(light_curve), dtype=[("size", float),("pop", float)])
	cell_curve['pop'] = light_curve['RATE'];

	# Retrieve the change points array (i.e. the edges of the Bayesian blocks)
	cp_array, num_cps, bestlogprob, lastcellstart = lcbayes(cell_curve, ncp_prior, nlag);
	print(cp_array)

	# Rebin light curve 
	cell_curve = rebinlc(light_curve, cp_array, num_cps)
	
	# Re-retrieve the change points
	cparray, num_cps, bestlogprob, lastcellstart = lcbayes(cell_curve, ncp_prior, nlag);
	
	# Compute the cumulative total number of counts, matched to the time array
	cumsum = lccumsum(light_curve);

	# Coalesce first or last blocks if they are outliers
	if ((num_cps > 2) and (coalescefrac > 0)):
		t0start = light_curve['TIME'][cparray[0]] - 0.5*timedel;
		t1start = light_curve['TIME'][cparray[1]] - 0.5*timedel;
		t0stop = t1start;
		t1stop = None;

		if (cparray[2] == nrows): 
			t1stop = light_curve['TIME'][nrows-1] + 0.5*timedel;
		else:
			t1stop = light_curve['TIME'][cparray[2]] - 0.5*timedel;

		# Remove the first changepoint
		if ((t0stop-t0start) < (t1stop-t1start)*coalescefrac):
			i = 1
			for i in range(num_cps-1):
				cparray[i] = cparray[i+1]
			num_cps-=1;


	if (num_cps > 2 and coalescefrac > 0):
		t0start = light_curve['TIME'][cparray[num_cps-3]] - 0.5*timedel;
		t1start = light_curve['TIME'][cparray[num_cps-2]] - 0.5*timedel;
		t0stop = t1start;
		t1stop = None;

		if (cparray[num_cps-1] == nrows):
			t1stop = light_curve['TIME'][nrows-1] + 0.5*timedel;
		else:
			t1stop  = light_curve['TIME'][cparray[num_cps-1]] - 0.5*timedel;

		# Remove the last changepoint
		if ((t1stop-t1start) < (t0stop-t0start)*coalescefrac):
			cparray[num_cps-2] = cparray[num_cps-1];
			num_cps-=1

	# Find bin edges. The last point is tricky, since it refers to the N+1'th cell. */
	bbgti = np.zeros(shape=num_cps-1, dtype=[("start",float),("stop",float)])
	for i in range(num_cps-1): 
		bbgti['start'][i] = light_curve['TIME'][cparray[i]] - 0.5*timedel;
		if (cparray[i+1] == nrows):
			bbgti["stop"][i]  = light_curve['TIME'][nrows-1] + 0.5*timedel;
		else:
			bbgti["stop"][i]  = light_curve['TIME'][cparray[i+1]] - 0.5*timedel;
	ngti = num_cps-1;

	# Estimate the burst duration, from the end of the first BB to the beginning of the last BB
	burst_tstart = bbgti["stop"][0];
	burst_tstop  = bbgti["start"][ngti-1];

	# Estimate the burst duration intervals
	duration = burst_tstop - burst_tstart

	# Find istart and istop
	istart, istop = burstspan(light_curve, burst_tstart, burst_tstop)

	# Estimate the fluence
	fluence = cumsum[istop] - cumsum[istart];

	print(burst_tstart)
	print(burst_tstop)

	return duration, fluence

def lcbayes(cell_curve, ncp_prior, nlag):
	"""
	Determine the Bayesian block change points, based on the 
	Poisson binned cost function
		
	Attributes:
	----------
	cell_curve : np.ndarray

	ncp_prior : float
		log parameter for number of changepoints
	"""

	ncells = len(cell_curve)
	cum_cell_curve = np.zeros(shape=ncells, dtype=[("size",float),("pop",float)])

	last_start = np.zeros(shape=ncells, dtype=int) # last cell start 
	merged = np.zeros(shape=ncells)
	best = np.zeros(shape=ncells) # best log probability 
	temp= 0;
	imaxer=0;

	istart = 0;
	ioldstart = 0;
	for i in range(ncells):
		# Approximation to the "nibble" algorithm
		if (nlag > 0):
			istart = i - nlag;
			if (istart < 0):
				istart = 0;

		# If we are nibbling, then we must shift the best[] array, so
		# that the normalized probability is unity before the starting
		# element.
		if ( (istart > 0) and (istart != ioldstart) ):
			for j in range(istart, i):
				best[j] -= best[istart-1];

		# Accumulate the parameters
		for j in range(istart, i):
			cum_cell_curve['size'][j] += cell_curve['size'][i];
			cum_cell_curve['pop'][j]  += cell_curve['pop'][i];
	

		cum_cell_curve['size'][i] = cell_curve['size'][i];
		cum_cell_curve['pop'][i] = cell_curve['pop'][i];

		# Compute the cost function for the cumulants */
		merged[istart:i+1-istart] = logprob_lc(cum_cell_curve[istart:i+1-istart], ncp_prior);

		# Where is the maximum probability in the joint best|mergedarrays? */
		imaxer = istart;
		best[i] = merged[istart];
		if (i > 0):
			for j in range(istart+1, i+1):
				temp = best[j-1]+merged[j];
				if (temp > best[i]):
					best[i] = temp;
					imaxer = j;


		# Record the new best position */
		last_start[i] = imaxer;

		# Keep track of the previous nibble starting point */
		ioldstart = istart;

	# Count number of change points */
	ncp = 2;
	index = last_start[ncells-1];
	while (index > 1):
		ncp +=1;
		index = last_start[index-1];

	# Create output array of change points
	cparray = np.zeros(shape=ncp, dtype=int)

	icp = ncp-1;
	cparray[icp] = ncells;
	icp -=1 # maybe this needs to go before cparray call
	index = last_start[ncells-1];
	while (index > 1):
		cparray[icp] = index;
		icp -=1
		index = last_start[index-1];
	cparray[0] = 0;



	return cparray, ncp, best, last_start


def logprob_lc(cell_curve, ncp_prior):
	"""
	Computes log posterior probability for binned data. See:  J.D. Scargle, 1998, ApJ, 504, 405.
	
	Attributes:
	----------
	cell_curve:

	ncp_prior : float
		log parameter for number of changepoints
	"""

	ncells = len(cell_curve)
	logprob = np.zeros(ncells)

	for i in range(ncells):
		logprob[i] = scipy.special.gammaln(cell_curve['pop'][i]+1) - (cell_curve['pop'][i]+1)*np.log(cell_curve['size'][i]);
		logprob[i] -= ncp_prior;

	return logprob

def burstspan(light_curve, tstart, tstop):
	"""
	Find start and stop times:
	1. start time is end of first BB, which is assumed to be the
	pre-burst background.
	2. stop is the start of the last BB, which is assumed to be the
	post-burst background.
	"""

	nrows = len(light_curve)

	for i in range(nrows):
		if (light_curve['TIME'][i] >= tstart):
			break;
	istart = i-1;
	
	if (i == 0):
		istart = 0;
	
	for j in range(istart+1,nrows):
		if (light_curve['TIME'][j] >= tstop):
			break;
	istop = j;
	if (i == nrows):
		istop = nrows-1;

	return istart, istop

def lccumsum(light_curve):
	"""
	Form the cumulative sum of the light curve, between two points.
	The data is assumed to be expressed in counts already. 
	Times are assumed to be center-bin.
	"""

	cumcounts = np.zeros(shape=len(light_curve));

	cumcounts[0] = light_curve['RATE'][0];
	for i in range(1, len(light_curve)):
		cumcounts[i] = cumcounts[i-1] + light_curve['RATE'][i];

	return cumcounts;



def rebinlc(light_curve, cp_array, num_cps):
	"""
	Rebin an existing light curve to a new one, given a set of change points
	"""

	timedel = light_curve['TIME'][1] - light_curve['TIME'][0]
	nc = num_cps-1;
	new_cell_curve = np.zeros(shape=nc,dtype=[("size",float),("pop",float)])
  
	for j in range(nc):
		dtot = 0;
		new_cell_curve["pop"][j] = 0;
		for i in range(cp_array[j], cp_array[j+1]):
			dtot += timedel # since we are using constant time bins, we can just use timedel here.
			new_cell_curve["pop"][j] += light_curve['RATE'][i];
		new_cell_curve['size'][j] = int(dtot/timedel);

	return new_cell_curve;
