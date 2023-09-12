"""
Author: Mike Moss
Contact: mikejmoss3@gmail.com	

Defines all the functions necessary to simulate an observation of a GRB using an input template, designated instrument response, and observing conditions 

"""

import numpy as np
from packages_simulations.GRB import GRB

from packages_util.cosmology import lum_dis
from packages_simulations.RSP import ResponseMatrix


def simulate_observation(template_grb, z, imx, imy, ndets, bgd, trials, resp_mat, multiproc=True, sim_triggers=False):
	"""
	Method to complete a simulation of a synthetic observation based on the input source frame GRB template and the desired observing conditions

	Attributes:
	template_grb = 		(GRB) GRB class object that holds the source frame information of the template GRB
	z =					(float) Redshift of synthetic GRB
	imx, imy = 			(float) The x and y position of the GRB on the detector plane
	ndets = 			(int) Number of detectors enabled during the synthetic observation 
	bdg = 				(float) Background level to be added to the synthetic light curve
	trials = 			(int) Number of trials to re-simulate this burst
	multiproc = 		(boolen) Multiprocessing on or off?
	sim_triggers = 		(boolean) Simulate triggers or not?
	"""

	# Initialize synth_GRB
	synth_GRB = GRB(grbname=template_grb.grbname,z=z,imx=imx,imy=imy)

	# Apply distance corrections to template GRB light curve to create synthetic GRB light cure
	synth_GRB.light_curve = template_grb.light_curve / 4. / np.pi / lum_dis(z)**2.

	# Apply observing condition corrections (e.g., NDETS)
	synth_GRB.light_curve /= ndets

	# Fold GRB through instrument response (RSP selected based on position on the detector plane)
	folded_spec = resp_mat.fold_spec(template_grb.spectrum)

	# Modulate the light curve by the folded spectrum normalization 

	# Add background to light curve 

	# Apply fluctuations 


	return synth_GRB

