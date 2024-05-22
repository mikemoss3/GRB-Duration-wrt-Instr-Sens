"""
Author: Mike Moss
Contact: mikejmoss3@gmail.com	

Test running sandbox and unit test runner

"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits 
import importlib
from util_packages.package_datatypes import dt_sim_res
from packages.package_bayesian_block import bayesian_t_blocks
from packages.class_SPECFUNC import PL, CPL
from packages.class_GRB import GRB
from packages.class_PLOTS import PLOTSIMRES, PLOTSAMPLE, PLOTGRB, PLOTS
from util_packages.package_det_ang_dependence import find_pcode, find_inc_ang, fraction_correction
from packages.package_many_simulations import many_simulations, make_param_list, make_ave_sim_res


grbs_names = np.array([
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
	"130427A_cut",
	"130603B",
	"130925A",
	"140506A",
	"160425A",
	"160804A",
	"161001A",
	"161219B",
	], dtype="U11")





for i in range(len(grbs_names)):
	grbp = importlib.import_module("data_files.grb_{}.info".format(grbs_names[i]), package=None) # Load GRB parameters
	
	template_grb = GRB(grbname = grbp.name, z=grbp.z)
	
	sim_results = np.load("data_files/grb_{}/grb_{}_redshift_sim-results.tmp.txt.npy".format(grbp.name, grbp.name, grbp.name))
	if template_grb.dt < 1:
		sim_results['FLUENCE'] *= 0.064
		sim_results['T100FLUENCE'] *= 0.064
		sim_results['1sPeakFlux'] *= 0.064
	np.save("data_files/grb_{}/grb_{}_redshift_sim-results.tmp.txt".format(grbp.name, grbp.name), sim_results)
	
