"""
Author: Mike Moss
Contact: mikejmoss3@gmail.com	

Defines the Bayesian block method to calculate the duration of a GRB from a supplied light curve

"""

import numpy as np
import scipy 
from astropy.stats import bayesian_blocks


def bayesian_t_blocks(grb,dur_per=90,ncp_prior=6):
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

	return 0, 0, 0

	return duration, t_start, phot_fluence


def custom_bblocks(grb):
	"""
	Bayesian block implementation based on heasoft tool battblocks
	"""
	# Bayesian blocks for pure Poisson counting data
	lc2cells(t, rate, dt, nrows, cellsizes, cellpops, ncells, timedel);

	# Retrieve the change points array (i.e. the edges of the Bayesian blocks)
	cparray = lcbayes(cellsizes, cellpops, ncells, ncp_prior, ncp, best, last_start, nlag);

	# Compute the cumulative total number of counts, matched to the time array
	cumsum = lccumsum(rate, nrows);

	# Coalesce first or last blocks if they are outliers
	if ((ncp > 2) and (coalescefrac > 0)):
		t0start = t[cparray[0]] - 0.5*dt[cparray[0]];
		t1start = t[cparray[1]] - 0.5*dt[cparray[1]];
		t0stop = t1start;
		t1stop = None;

		if (cparray[2] == nrows): 
			t1stop = t[nrows-1] + 0.5*dt[nrows-1];
		else:
			t1stop = t[cparray[2]] - 0.5*dt[cparray[2]];

		# Remove the first changepoint
		if ((t0stop-t0start) < (t1stop-t1start)*coalescefrac):
			i = 1
			for i in range(ncp-1):
				cparray[i] = cparray[i+1]
			ncp-=1;


	if (ncp > 2 and coalescefrac > 0):
		t0start = t[cparray[ncp-3]] - 0.5*dt[cparray[ncp-3]];
		t1start = t[cparray[ncp-2]] - 0.5*dt[cparray[ncp-2]];
		t0stop = t1start;
		t1stop = None;

		if (cparray[ncp-1] == nrows):
			t1stop = t[nrows-1] + 0.5*dt[nrows-1];
		else:
			t1stop  = t[cparray[ncp-1]] - 0.5*dt[cparray[ncp-1]];

		# Remove the last changepoint
		if ((t1stop-t1start) < (t0stop-t0start)*coalescefrac):
			cparray[ncp-2] = cparray[ncp-1];
			ncp-=1

	# Find bin edges. The last point is tricky, since it refers to the N+1'th cell. */
	for i in range(ncp-1): 
		bbgti.start[i] = t[cparray[i]] - 0.5*dt[cparray[i]];
		if (cparray[i+1] == nrows):
			bbgti.stop[i]  = t[nrows-1] + 0.5*dt[nrows-1];
		else:
			bbgti.stop[i]  = t[cparray[i+1]] - 0.5*dt[cparray[i+1]];
	bbgti.ngti = ncp-1;

	# Estimate the burst duration, from the end of the first BB to the beginning of the last BB
	burst_tstart = bbgti.stop[0];
	burst_tstop  = bbgti.start[bbgti.ngti-1];

	# Find istart and istop
	burstspan(t, dt, nrows, burst_tstart, burst_tstop, istart, istop)

	# Estimate the fluence
	fluence = cumsum[istop] - cumsum[istart];


def burstspan(t):
	"""
	Find start and stop times:
	1. start time is end of first BB, which is assumed to be the
	pre-burst background.
	2. stop is the start of the last BB, which is assumed to be the
	post-burst background.
	"""
	istart, istop = None, None;
	i = 0;
	j = 0;

	for i in range(nrows):
		if (t[i] >= tstart):
			break;
	istart = i-1;
	
	if (i == 0):
		istart = 0;
	
	for j in range(istart+1,nrows):
		if (t[j] >= tstop):
			break;
	istop = j;
	if (i == nrows):
		istop = nrows-1;

	return istart, istop

def lccumsum(counts, ntimes):
	"""
	Form the cumulative sum of the light curve, between two points.
	The data is assumed to be expressed in counts already. 
	Times are assumed to be center-bin.
	"""

	cumcounts = np.zeros(shape=ntimes);

	cumcounts[0] = counts[0];
	for i in range(1, ntimes):
		cumcounts[i] = cumcounts[i-1] + counts[i];

	return cumcounts;

def lcbayes(cellsizes, cellpops, ncells, 
	ncp_prior, ncparray, bestlogprob, lastcellstart, nlag):
	"""
	Determine the Bayesian block change points, based on the 
	Poisson binned cost function
		
	Attributes:
	----------
	cellsizes : np.ndarray
		widths of cells, in units of timedel
	cellpops : np.ndarray
		number of events in each cell
	ncells : float 
		number of cells
	ncp_prior : float
		log parameter for number of changepoints
	ncparray :

	bestlogprob :

	lastcellstart :

	nlag :

	"""

	cumsizes = np.zeros(shape=ncells)
	cumpops = np.zeros(shape=ncells)
	last_start = np.zeros(shape=ncells)
	cparray = 0;
	merged = np.zeros(shape=ncells)
	best = np.zeros(shape=ncells)
	temp= 0;
	imaxer, ncp, index, icp = 0,0,0,0;
	istart = 0
	ioldstart = 0;

	if (bestlogprob):
		bestlogprob = 0;
	if (lastcellstart):
		lastcellstart = 0;

	for i in range(ncells):
		cumsizes[i] = 0;
		cumpops[i] = 0;

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

			# Accumulate the parameters */
			for j in range(istart, i):
				cumsizes[j] += cellsizes[i];
				cumpops[j]  += cellpops[i];
		
			cumsizes[i] = cellsizes[i];
			cumpops[i] = cellpops[i];

			# Compute the cost function for the cumulants */
			logprob_lc(cumpops+istart, cumsizes+istart, i+1-istart, merged+istart, ncp_prior);

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

	# Create output array of change points */
	cparray = np.zeros(shape=ncp);
	if (cparray == 0):
		ncp = 0;

	icp = ncp-1;
	cparray[icp] = ncells;
	icp -=1 # maybe this needs to go before cparray call
	index = last_start[ncells-1];
	while (index > 1):
		cparray[icp] = index;
		icp -=1
		index = last_start[index-1];
	cparray[0] = 0;

	ncparray = ncp;
	return cparray;


def logprob_lc(cellpops, cellsizes, ncells, ncp_prior):
	"""
	Computes log posterior probability for binned data. See:  J.D. Scargle, 1998, ApJ, 504, 405.
	
	Attributes:
	----------
	cellpops : np.ndarray
		number of events in each cell
	cellsizes : np.ndarray
		widths of cells, in units of timedel
	ncells : float 
		number of cells
	ncp_prior : float
		log parameter for number of changepoints
	"""

	logprob = np.zeros(ncells)

	for i in range(ncells):
		logprob[i] = scipy.special.lgamma(cellpops[i]+1) - (cellpops[i]+1)*np.log(cellsizes[i]);
		logprob[i] -= ncp_prior;

	return logprob

def lc2cells(t, counts, dt, ntimes, cellsizes, cellpops, ncells, timedel):
	"""
	Convert a light curve to cells, by direct transcription

	Attributes:
	----------
	t : np.ndarray
		array of light curve times
	counts : np.ndarray
		array of light curve counts
	dt : np.ndarry
		array of time bin sizes
	cellsizes : np.ndarray
		widths of cells, in units of timedel
	cellpops : np.ndarray
		number of events in each cell
	ncells : float 
		number of cells
	timedel : float 
		time bin size
	"""
	cpops = 0
	csize = 0

	if ((t == 0) or (counts == 0) or (dt == 0) or 
		(ntimes <= 0) or (cellsizes == 0) or 
		(cellpops == 0)):
		return 0;

	nc = ntimes;

	if nc == 0:
		return 0;

	for i in range(nc):
		csize[i] = dt[i] / timedel;
		cpops[i] = counts[i];

	cellsizes = csize;
	cellpops = cpops;
	ncells = nc;

	return 0;