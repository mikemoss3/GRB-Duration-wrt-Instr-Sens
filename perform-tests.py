"""
Author: Mike Moss
Contact: mikejmoss3@gmail.com	

Test running sandbox and unit test runner

"""

import numpy as np
from astropy.io import fits 
import importlib
from util_packages.package_datatypes import dt_sim_res
from packages.package_bayesian_block import bayesian_t_blocks
from packages.class_GRB import GRB
from packages.class_PLOTS import PLOTSAMPLE

low_z_grbs_names = np.array([
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

high_z_grbs_names = np.array([
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



low_z_results = np.zeros(shape=0, dtype=dt_sim_res)
low_z_results_unshift = np.zeros(shape=0, dtype=dt_sim_res)
for i in range(len(low_z_grbs_names)):
	grbp = importlib.import_module("data_files.grb_{}.info".format(low_z_grbs_names[i]), package=None) # Load GRB parameters
	sim_results = np.load("data_files/grb_{}/grb_{}_redshift_sim-results.tmp.txt.npy".format(grbp.name, grbp.name))
	
	low_z_results = np.append(low_z_results, sim_results[(sim_results["z"]>3) & (sim_results["z"]<9)])

	low_z_results_unshift = np.append(low_z_results_unshift, sim_results[sim_results["z"]<1])


obs_lobs_z_results = np.zeros(shape=len(low_z_grbs_names), dtype=[("DURATION",float),("z",float)])
for i in range(len(low_z_grbs_names)):
	grbp = importlib.import_module("data_files.grb_{}.info".format(low_z_grbs_names[i]), package=None) # Load GRB parameters
	obs_lobs_z_results[i]["DURATION"] = grbp.t_true
	obs_lobs_z_results[i]["z"] = grbp.z

high_z_results = np.zeros(shape=len(high_z_grbs_names), dtype=[("DURATION",float),("z",float)])
for i in range(len(high_z_grbs_names)):
	grbp = importlib.import_module("data_files.grb_{}.info".format(high_z_grbs_names[i]), package=None) # Load GRB parameters
	high_z_results[i]["DURATION"] = grbp.t_true
	high_z_results[i]["z"] = grbp.z


sim_res_plots = PLOTSAMPLE(data_table=low_z_results)
sim_res_plots.cumulative_durations()
sim_res_plots.show()