"""
Author: Mike Moss
Contact: mikejmoss3@gmail.com	

Defines the functions that generate the duration and fluence distributions of synthetic GRBs
"""

import numpy as np
import multiprocessing as mp
from packages.package_simulations import simulate_observation
from packages.package_bayesian_block import bayesian_t_blocks


def many_simulations(template_grb, resp_mat, z_arr, imx_arr, imy_arr, ndets_arr, trials, dur_per = 90, multiproc=True, num_cores = 4, sim_triggers=False,ndet_max=32768,bgd_rate_per_det=0.3):
	"""
	Method to perform multiple simulations for each combination of input parameters 

	Attributes:
	----------
	template_grb : GRB
		GRB object used as a template light curve and spectrum
	resp_mat : RSP
		Response matrix to convolve the template spectrum with 
	z_arr : np.ndarry
		Array of redshifts to simulate the GRB at
	imx_arr, imy_arr : np.ndarry, np.ndarry
		Array of (imx,imy) points i.e., where the simualted sources will be located on the detector 
	ndets_arr : np.ndarry
		Arry of values to use for the number of enabled detectors during the observation simulations
	trials : int
		Number of trials for each parameter combination
	dur_per : float
		Duration percentage to find using Bayesian blocks, i.e., dur_pur = 90 returns the T_90 of the burst
	multiproc : bool
		Whether to multiprocessing or not
	num_cores : int
		Number of cores to use for multiprocessing (this must be <= the number of cpu's the computer has)
	sim_triggers : bool
		Whether or not to simulate the Swift/BAT trigger algorithms or not
	ndet_max : int
		Maximum number of detectors on the detector plane (for Swift/BAT ndet_max = 32,768)
	bgd_rate_per_det : float
		Background rate per detector (for Swift/BAT bgd_rate_per_det is close to 0.3 cnts / s / det across mission lifetime)
	"""

	# Make a list of all parameter combinations	
	param_list = np.array(np.meshgrid(z_arr, imx_arr, imy_arr,ndets_arr)).T.reshape(-1,4)

	sim_results = np.zeros(shape=int(len(param_list)*trials),dtype=[("DURATION",float),("TSTART",float),("FLUENCE",float),("z",float),("imx",float),("imy",float),("ndets",float)])
	sim_result_ind = 0

	# Simulate an observation for each parameter combination
	for i in range(len(param_list)):

		if multiproc is False:
			# Run the simulations without multiprocessing 
			for j in range(trials):

				sim_results[["z", "imx", "imy", "ndets"]][sim_result_ind] = (param_list[i][0], param_list[i][1], param_list[i][2], param_list[i][3])

				synth_GRB = simulate_observation(template_grb=template_grb,z_p=param_list[i][0],imx=param_list[i][1],imy=param_list[i][2],ndets=param_list[i][3],resp_mat=resp_mat,sim_triggers=sim_triggers,ndet_max=ndet_max,bgd_rate_per_det=bgd_rate_per_det)
				sim_results[["DURATION", "TSTART", "FLUENCE"]][sim_result_ind] = bayesian_t_blocks(synth_GRB, dur_per=dur_per) # Find the Duration and the fluence 

				# Increase simulation index
				sim_result_ind +=1
		else:
			# Run the simulations with multiprocessing
			 
			# Load in a number of pools to run the code.
			with mp.Pool(num_cores) as pool:
				synth_GRBs = pool.map(simulate_observation, [template_grb, param_list[i][0], param_list[i][1], param_list[i][2], param_list[i][3], resp_mat,sim_triggers, ndet_max, bgd_rate_per_det] )

			# Add the new results to the list of sim results
			for k in range(num_cores):

				sim_results[["z","imx","imy","ndets"]][sim_result_ind] = (param_list[i][0], param_list[i][1], param_list[i][2], param_list[i][3])

				sim_results[["DURATION", "TSTART", "FLUENCE"]][sim_result_ind] = bayesian_t_blocks(t=synth_GRBs[k], dur_per=dur_per)
				sim_result_ind += k
						
			# Increment the index tacking value by the number of trials just simulated
			# sim_result_ind += trials

	return sim_results

