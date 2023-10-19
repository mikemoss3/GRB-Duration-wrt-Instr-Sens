"""
Author: Mike Moss
Contact: mikejmoss3@gmail.com	

Defines functions to obtain the duration and fluence measurements for many synthetic GRBs
"""

import numpy as np
import multiprocessing as mp
from packages.class_GRB import GRB
from packages.class_RSP import ResponseMatrix
from packages.package_simulations import simulate_observation
from packages.package_bayesian_block import bayesian_t_blocks
from util_packages.package_datatypes import dt_sim_res

import matplotlib.pyplot as plt

def many_simulations(template_grb, param_list, trials, dur_per = 90, 
	multiproc=True, num_cores = 4, sim_triggers=False, ndet_max=32768, bgd_rate_per_det=0.3, 
	out_file_name = None, ret_ave = False, keep_synth_grbs=False):
	"""
	Method to perform multiple simulations for each combination of input parameters 

	Attributes:
	----------
	template_grb : GRB
		GRB object used as a template light curve and spectrum
	z_arr : np.ndarray
		Array of redshifts to simulate the GRB at
	imx_arr, imy_arr : np.ndarry, np.ndarray
		Array of (imx,imy) points i.e., where the simualted sources will be located on the detector 
	ndets_arr : np.ndarray
		Array of values to use for the number of enabled detectors during the observation simulations
	trials : int
		Number of trials for each parameter combination
	dur_per : float
		Duration percentage to find using Bayesian blocks, i.e., dur_pur = 90 returns the T_90 of the burst
	multiproc : boolean
		Whether to multiprocessing or not
	num_cores : int
		Number of cores to use for multiprocessing (this must be <= the number of cpu's the computer has)
	sim_triggers : boolean
		Whether or not to simulate the Swift/BAT trigger algorithms or not
	ndet_max : int
		Maximum number of detectors on the detector plane (for Swift/BAT ndet_max = 32,768)
	bgd_rate_per_det : float
		Background rate per detector (for Swift/BAT bgd_rate_per_det is close to 0.3 cnts / s / det across mission lifetime)
	out_file_name : string
		File name of .txt file to write the simulation result out to. 
	"""

	# Make a list to hold the simulation results
	sim_results = np.zeros(shape=int(len(param_list)*trials),dtype=dt_sim_res)
	sim_result_ind = 0

	# Make a Response Matrix object
	resp_mat = ResponseMatrix()

	if keep_synth_grbs is True:
		synth_grb_arr = np.zeros(shape=len(param_list),dtype=GRB)

	# Simulate an observation for each parameter combination
	for i in range(len(param_list)):

		if multiproc is False:
			# Run the simulations without multiprocessing 
			for j in range(trials):

				sim_results[["z", "imx", "imy", "ndets"]][sim_result_ind] = (param_list[i][0], param_list[i][1], param_list[i][2], param_list[i][3])

				# Load Swift BAT response based on position on detector plane  
				resp_mat.load_SwiftBAT_resp(sim_results[sim_result_ind]["imx"],sim_results[sim_result_ind]["imy"])


				synth_GRB = simulate_observation(template_grb=template_grb,z_p=param_list[i][0],imx=param_list[i][1],imy=param_list[i][2],ndets=param_list[i][3],resp_mat=resp_mat,sim_triggers=sim_triggers,ndet_max=ndet_max,bgd_rate_per_det=bgd_rate_per_det)
				sim_results[["DURATION", "TSTART", "FLUENCE"]][sim_result_ind] = bayesian_t_blocks(synth_GRB, dur_per=dur_per) # Find the Duration and the fluence 

				# Increase simulation index
				sim_result_ind +=1
			if keep_synth_grbs is True:
				synth_grb_arr[i] = synth_GRB
		else:
			# Run the simulations with multiprocessing
			 
			# Load in a number of pools to run the code.
			with mp.Pool(num_cores) as pool:
				synth_GRBs = pool.starmap(simulate_observation, [(template_grb, param_list[i][1], param_list[i][2], param_list[i][3], None, param_list[i][0], sim_triggers, ndet_max, t) for t in range(trials)])

			# Add the new results to the list of sim results
			for k in range(trials):

				sim_results[["z","imx","imy","ndets"]][sim_result_ind] = (param_list[i][0], param_list[i][1], param_list[i][2], param_list[i][3])

				sim_results[["DURATION", "TSTART", "FLUENCE"]][sim_result_ind] = bayesian_t_blocks(synth_GRBs[k], dur_per=dur_per) # Find the Duration and the fluence 
				sim_result_ind += 1


			if keep_synth_grbs is True:
				synth_grb_arr[i] = synth_GRBs[0]

						
			# Increment the index tacking value by the number of trials just simulated
			# sim_result_ind += trials

	if out_file_name is not None:
		np.savetxt(out_file_name,sim_results)

	if ret_ave is True:
		sim_results = make_ave_sim_res(sim_results)

	if keep_synth_grbs is True:
		return sim_results, synth_grb_arr
	else:
		return sim_results

def make_param_list(z_arr, imx_arr, imy_arr, ndets_arr):
	"""
	Method to make a list of all parameter combinations from the given parameter values.

	Attributes:
	----------
	z_arr : np.ndarray
		Array of redshifts to simulate the GRB at
	imx_arr, imy_arr : np.ndarry, np.ndarray
		Array of (imx,imy) points i.e., where the simualted sources will be located on the detector 
	ndets_arr : np.ndarray
		Array of values to use for the number of enabled detectors during the observation simulations
	"""

	# Make a list of all parameter combinations	
	param_list = np.array(np.meshgrid(z_arr, imx_arr, imy_arr,ndets_arr)).T.reshape(-1,4)

	return param_list


def make_ave_sim_res(sim_results):
	"""
	Method to make an average sim_results array for each unique parameter combination

	Attributes:
	----------
	sim_results : np.ndarray
		sim_results array 
	"""

	unique_rows = np.unique(sim_results[["z","imx","imy","ndets"]])

	ave_sim_results = np.zeros(shape=len(unique_rows),dtype=dt_sim_res)

	ave_sim_results[["z","imx","imy","ndets"]] = unique_rows

	for i in range(len(unique_rows)):
		ave_sim_results["DURATION"][i] = np.sum(sim_results["DURATION"][ sim_results[["z","imx","imy","ndets"]] == unique_rows[i]])/len(sim_results[ sim_results[["z","imx","imy","ndets"]] == unique_rows[i]])
		ave_sim_results["TSTART"][i] = np.sum(sim_results["TSTART"][ sim_results[["z","imx","imy","ndets"]] == unique_rows[i]])/len(sim_results[ sim_results[["z","imx","imy","ndets"]] == unique_rows[i]])
		ave_sim_results["FLUENCE"][i] = np.sum(sim_results["FLUENCE"][ sim_results[["z","imx","imy","ndets"]] == unique_rows[i]])/len(sim_results[ sim_results[["z","imx","imy","ndets"]] == unique_rows[i]])

	return ave_sim_results