"""
Author: Mike Moss
Contact: mikejmoss3@gmail.com	

Defines the functions that generate the duration and fluence distributions of synthetic GRBs
"""

import numpy as np
import multiprocessing as mp
from packages.package_simulations import simulate_observation
from packages.package_bayesian_block import bayesian_t_blocks


def many_simulations(template_grb, resp_mat, z_arr, imx_arr, imy_arr, ndets_arr, trials, multiproc=True, sim_triggers=False,ndet_max=32768,bgd_rate_per_det=0.3,num_cores = 4):
	"""
	Method to perform multiple simulations for each combination of input parameters 

	Attributes:

	num_cores = Number of cores to use for multiprocessing (this must be <= the number of cpu's the computer has)
	"""


	# Make a list of all parameter combinations	
	param_list = np.array(np.meshgrid(z_arr, imx_arr, imy_arr,ndets_arr)).T.reshape(-1,4)

	sim_results = np.zeros(shape=int(len(param_list)*trials),dtype=[("T90",float),("FLUENCE",float)])
	sim_result_ind = 0

	# Simulate an observation for each parameter combination
	for i in range(len(param_list)):

		if multiproc is False:
			# Run the simulations without multiprocessing 
			for j in range(trials):

				synth_GRB = simulate_observation(template_grb=template_grb,z=param_list[i][0],imx=param_list[i][1],imy=param_list[i][2],ndets=param_list[i][3],resp_mat=resp_mat,sim_triggers=sim_triggers,ndet_max=ndet_max,bgd_rate_per_det=bgd_rate_per_det)
				sim_results[sim_result_ind] = bayesian_t_blocks(synth_GRB) # Find the T90 and the fluence 

				# Increase simulation index
				sim_result_ind +=1
		else:
			# Run the simulations with multiprocessing
			 
			# Load in a number of pools to run the code.
			with mp.Pool(num_cores) as pool:
				synth_GRBs = pool.map(simulate_observation, [template_grb,param_list[i][0],param_list[i][1],param_list[i][2],param_list[i][3],resp_mat,sim_triggers,ndet_max,bgd_rate_per_det] )

			# Add the new results to the list of sim results
			for k in range(num_cores):
				sim_results[sim_result_ind] = bayesian_t_blocks(t=synth_GRBs[k].light_curve['TIME'],x=synth_GRBs[k].light_curve['RATE'],sigma=synth_GRBs[k].light_curve['UNC'],fitness="measures")
				sim_result_ind += k
						
			# Increment the index tacking value by the number of trials just simulated
			# sim_result_ind += trials

	return sim_results

