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
	# "050416A",
	# "050525A",
	# "060614",
	# "060912A",
	"061021",
	# "080430",
	# "080916A",
	# "081007",
	# "090424",
	# "091018",
	# "091127",
	# "100621A",
	# "100625A",
	# "100816A",
	# "101219A",
	# "110715A",
	# "111228A",
	# "120311A",
	# "130427A",
	# "130427A_cut",
	# "130603B",
	# "130925A",
	# "140506A",
	# "160425A",
	# "160804A",
	# "161001A",
	# "161219B",
	], dtype="U11")



grbp = importlib.import_module("data_files.grb_{}.info".format(grbs_names[0]), package=None) # Load GRB parameters

template_grb = GRB(grbname = grbp.name, z=grbp.z)
template_grb.load_light_curve(grbp.fn, rm_trigtime=True, det_area=0.16)
template_grb.cut_light_curve(tmin=grbp.t_cut_min, tmax=grbp.t_cut_max, buffer=grbp.t_buffer)
template_grb.load_specfunc(CPL(alpha= grbp.alpha, ep=grbp.ep, norm=grbp.norm, enorm=50))

z_arr = np.array([grbp.z])
imx_arr = np.array([0.])
imy_arr = np.array([0.])
ndets_arr = np.array([30000])

param_list = make_param_list(z_arr,imx_arr,imy_arr,ndets_arr)
trials = 1

sim_results, synth_grbs = many_simulations(template_grb, param_list, trials, multiproc=False, keep_synth_grbs=True, verbose=True)

plot = PLOTGRB()
ax = plt.figure().gca()
plot.plot_light_curves(grbs=template_grb, ax=ax)
plot.plot_light_curves(grbs=synth_grbs, ax=ax)

plot.show()