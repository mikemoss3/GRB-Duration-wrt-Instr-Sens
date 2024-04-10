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

# Grab obs high-z bursts
obs_high_z = np.zeros(shape=len(obs_high_z_grbs),dtype=[("DURATION",float)])
for i in range(len(obs_high_z_grbs)):
	grbp = importlib.import_module("data_files.grb_{}.info".format(obs_high_z_grbs[i]), package=None) # Load GRB parameters
	obs_high_z['DURATION'][i] = grbp.t_true

res_lowobs_lowsim = ks_2samp(obs_low_z['DURATION'], sim_low_z['DURATION'])
print(res_lowobs_lowsim)
res_lowsim_highsim = ks_2samp(sim_low_z['DURATION'], sim_high_z['DURATION'])
print(res_lowsim_highsim)
res_lowsim_highobs = ks_2samp(sim_low_z['DURATION'], obs_high_z['DURATION'])
print(res_lowsim_highobs)
res_highobs_highsim = ks_2samp(obs_high_z['DURATION'], sim_high_z['DURATION'])
print(res_highobs_highsim)
res_highobs_lowobs = ks_2samp(obs_high_z['DURATION'], obs_low_z['DURATION'])
print(res_highobs_lowobs)