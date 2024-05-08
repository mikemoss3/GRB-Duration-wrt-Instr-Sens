import importlib
import numpy as np
from scipy.stats import ks_2samp, kstest

from util_packages.package_datatypes import dt_sim_res

obs_low_z_grbs = np.array([
	"050416A",
	"050525A",
	"060614",
	"060912A",
	"061021",
	"080430",
	"080916A",
	"081007",
	"090424",
	"091018",
	"091127",
	"100621A",
	"100625A",
	"100816A",
	"101219A",
	"110715A",
	"111228A",
	"120311A",
	"130427A",
	"130603B",
	"130925A",
	"140506A",
	"160425A",
	"160804A",
	"161001A",
	"161219B",
	], dtype="U10")

obs_high_z_grbs = np.array([
	"060206", 
	"060210", 
	"060306", 
	"060927", 
	"080607", 
	"090715B", 
	"111008A",
	"120712A",
	"130408A",
	"130606A",
	"170202A",
	], dtype="U10")

def calc_ks_test():
	# Grab obs low-z bursts
	obs_low_z = np.zeros(shape=len(obs_low_z_grbs),dtype=[("DURATION",float)])
	for i in range(len(obs_low_z_grbs)):
		grbp = importlib.import_module("data_files.grb_{}.info".format(obs_low_z_grbs[i]), package=None) # Load GRB parameters
		obs_low_z['DURATION'][i] = grbp.t_true

	# Grab simualted low-z bursts
	sim_low_z = np.zeros(shape=0, dtype=dt_sim_res)
	for i in range(len(obs_low_z_grbs)):
		grbp = importlib.import_module("data_files.grb_{}.info".format(obs_low_z_grbs[i]), package=None) # Load GRB parameters
		sim_results = np.load("data_files/results_final/grb_{}_redshift_sim-results.txt.npy".format(grbp.name, grbp.name))
		sim_low_z = np.append(sim_low_z, sim_results[sim_results["z"]<1])
		# sim_low_z = np.append(sim_low_z, sim_results)

	# Grab simualted high-z bursts
	sim_high_z = np.zeros(shape=0, dtype=dt_sim_res)
	for i in range(len(obs_low_z_grbs)):
		grbp = importlib.import_module("data_files.grb_{}.info".format(obs_low_z_grbs[i]), package=None) # Load GRB parameters
		sim_results = np.load("data_files/results_final/grb_{}_redshift_sim-results.txt.npy".format(grbp.name, grbp.name))
		sim_high_z = np.append(sim_high_z, sim_results[(sim_results["z"]>3) & (sim_results["z"]<9)])


	obs_high_z_peak_fluxes = np.array([0.3696, 0.167, 0.7110, 0.356, 2.927, 0.487, 0.808, 0.298, 0.740, 0.3409, 0.6056,])

	# Grab obs high-z bursts
	obs_high_z = np.zeros(shape=len(obs_high_z_grbs),dtype=[("DURATION",float), ("1sPeakFlux",float)])
	for i in range(len(obs_high_z_grbs)):
		grbp = importlib.import_module("data_files.grb_{}.info".format(obs_high_z_grbs[i]), package=None) # Load GRB parameters
		obs_high_z['DURATION'][i] = grbp.t_true
		obs_high_z['1sPeakFlux'][i] = obs_high_z_peak_fluxes[i]

	print("Durations:")
	res_lowobs_lowsim = ks_2samp(obs_low_z['DURATION'], sim_low_z['DURATION']).pvalue
	print("low-z obs + low-z sim: ",res_lowobs_lowsim)
	res_lowsim_highsim = ks_2samp(sim_low_z['DURATION'], sim_high_z['DURATION']).pvalue
	print("low-z sim + high-z sim: ",res_lowsim_highsim)
	res_highobs_lowobs = ks_2samp(obs_high_z['DURATION'], obs_low_z['DURATION']).pvalue
	print("high-z obs + low-z obs: ",res_highobs_lowobs)
	res_lowsim_highobs = ks_2samp(sim_low_z['DURATION'], obs_high_z['DURATION']).pvalue
	print("low-z sim + high-z obs: ",res_lowsim_highobs)
	res_highobs_highsim = ks_2samp(obs_high_z['DURATION'], sim_high_z['DURATION']).pvalue
	print("high-z obs + high-z sim: ",res_highobs_highsim)
	# print("\nFluences:")
	# res_lowobs_lowsim = ks_2samp(obs_low_z['1sPeakFlux'], sim_low_z['1sPeakFlux']).pvalue
	# print("low-z obs + low-z sim: ",res_lowobs_lowsim)
	# res_lowsim_highsim = ks_2samp(sim_low_z['1sPeakFlux'], sim_high_z['1sPeakFlux']).pvalue
	# print("low-z sim + high-z sim: ",res_lowsim_highsim)
	# res_highobs_lowobs = ks_2samp(obs_high_z['1sPeakFlux'], obs_low_z['1sPeakFlux']).pvalue
	# print("high-z obs + low-z obs: ",res_highobs_lowobs)
	# res_lowsim_highobs = ks_2samp(sim_low_z['1sPeakFlux'], obs_high_z['1sPeakFlux']).pvalue
	# print("low-z sim + high-z obs: ",res_lowsim_highobs)
	# res_highobs_highsim = ks_2samp(obs_high_z['1sPeakFlux'], sim_high_z['1sPeakFlux']).pvalue
	# print("high-z obs + high-z sim: ",res_highobs_highsim)

def calc_z_dist():
	sim_high_z = np.zeros(shape=0, dtype=dt_sim_res)
	z_array = np.zeros(shape=17, dtype=[("z_edges",float),("numGRBs",int)])
	z_array['z_edges'] = np.linspace(3, 11, num=len(z_array))
	print("z_min - z_max, GRB Remaining/Total")
	for i in range(len(obs_low_z_grbs)):
		grbp = importlib.import_module("data_files.grb_{}.info".format(obs_low_z_grbs[i]), package=None) # Load GRB parameters
		sim_results = np.load("data_files/results_final/grb_{}_redshift_sim-results.txt.npy".format(grbp.name, grbp.name))
		unique_z = np.unique(sim_results['z'])
		for j in range(len(z_array)-1):
			if any(unique_z > z_array['z_edges'][j]):
				z_array['numGRBs'][j] += 1
	

	for j in range(len(z_array)-1):
		print("{:.2f} - {:.2f}, {}/{} ".format(z_array['z_edges'][j], z_array['z_edges'][j+1], z_array['numGRBs'][j], len(obs_low_z_grbs)))

if __name__ == "__main__":
	calc_z_dist()
