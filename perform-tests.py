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
	# "130427A",
	"130603B",
	"130925A",
	"140506A",
	"160425A",
	"160804A",
	"161001A",
	"161219B",
	], dtype="U10")



for i in range(len(obs_low_z_grbs)):
	grbp = importlib.import_module("data_files.grb_{}.info".format(obs_low_z_grbs[i]), package=None) # Load GRB parameters
	sim_results = np.load("data_files/results_final/grb_{}_redshift_sim-results.txt.npy".format(grbp.name, grbp.name))
	if grbp.fn[-7:] == "64ms.lc":
		sim_results['1sPeakFlux']*=0.064

	np.save("data_files/results_final/grb_{}_redshift_sim-results.txt.npy".format(grbp.name), sim_results)
	