"""
Author: Mike Moss
Contact: mikejmoss3@gmail.com	

Defines the Bayesian block method to calculate the duration of a GRB from a supplied light curve

"""

import numpy as np
import astropy

def bayesian_blocks(grb,dur_per=90):
	"""
	Method to extract the duration of a GRB from a supplied light curve using a Bayesian block method. 

	Attributes:
	grb = 		(GRB) 		A grb object 
	dur_per =	(float) 	Percentage of the fluence to calculate the duration for (i.e., T90 corresponds to dur_per = 90)
	"""

	# Bin the light curve 
	bin_edges = astropy.stats.bayesian_block(t=synth_GRB.light_curve['TIME'],x=synth_GRB.light_curve['RATE'],sigma=synth_GRB.light_curve['UNC'],fitness="measures") # Find the T90 and the fluence 

	# Extract the duration
	# perhaps the bin_edges from above can be used to find the boxes that exceed some rate threshold?

	# return duration, fluence