"""
Author: Mike Moss
Contact: mikejmoss3@gmail.com	

Sandbox to perform all simulation, analysis, and plotting commands

"""

import importlib
import numpy as np 
import matplotlib.pyplot as plt
import time
from datetime import date

from packages.class_GRB import GRB
from packages.class_PLOTS import PLOTSIMRES, PLOTSAMPLE
from packages.class_SPECFUNC import PL, CPL
from packages.package_many_simulations import many_simulations, make_param_list, make_ave_sim_res
from util_packages.package_datatypes import dt_sim_res

def make_template_grb(grbp):

	## Make GRB object
	template_grb = GRB(grbname = grbp.name, z=grbp.z)
	# Load light curve 
	template_grb.load_light_curve(grbp.fn, rm_trigtime=True, det_area=0.16)
	template_grb.cut_light_curve(tmin=grbp.t_cut_min, tmax=grbp.t_cut_max)
	# Load spectrum
	template_grb.load_specfunc(CPL(alpha= grbp.alpha, ep=grbp.ep, norm=grbp.norm, enorm=50))

	return template_grb

def make_param_space(grbp):

	## Make parameter space  
	# z_arr = np.array([grbp.z])
	z_arr = np.linspace(grbp.z, grbp.zmax, num=50)
	imx_arr = np.array([0.])
	imy_arr = np.array([0.])
	# imx_arr = np.linspace(-1.75,1.75,70)
	# imy_arr = np.linspace(-0.875,0.875,70)
	ndets_arr = np.array([30000])

	param_list = make_param_list(z_arr,imx_arr,imy_arr,ndets_arr)
	return param_list

def main(name, template_grb, param_list, trials):

	sim_results = many_simulations(template_grb, param_list, trials, multiproc=False, keep_synth_grbs=False, verbose=True)
	ave_sim_results = make_ave_sim_res(sim_results)

	np.save("data_files/grb_{}/grb_{}_redshift_sim-results.tmp.txt".format(name, name), sim_results)
	np.save("data_files/grb_{}/grb_{}_redshift_ave-sim-results.tmp.txt".format(name, name), ave_sim_results)

	# np.save("data-files/grb-{}/grb_{}_detector_sim-results.tmp.txt".format(name, name), sim_results)
	# np.save("data-files/grb-{}/grb_{}_detector_ave-sim-results.tmp.txt".format(name, name), ave_sim_results)
	
	# sim_results, synth_grbs = many_simulations(template_grb, param_list, trials, multiproc=False, keep_synth_grbs=True, verbose=True)
	# return synth_grbs

def plot(name, t_true):
	sim_results = np.load("data_files/grb_{}/grb_{}_redshift_sim-results.tmp.txt.npy".format(name, name))
	plot = PLOTSIMRES() # Plot simulation results
	plot.redshift_evo(sim_results, t_true=t_true, log=False)
	
	# today = date.today()
	# plt.savefig("data_files/figs/z-evo-plots/grb-{}-redshift-evo.png".format(name), dpi=400)
	
	# plt.show()


if __name__ == "__main__":

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
		"130603B",
		"130925A",
		"140506A",
		"160425A",
		"160804A",
		"161001A",
		"161219B",
		], dtype="U10")

	"""
	for i in range(len(grbs_names)):
		grbp = importlib.import_module("data_files.grb_{}.info".format(grbs_names[i]), package=None) # Load GRB parameters

		template_grb = make_template_grb(grbp) # Create template GRB

		param_list = make_param_space(grbp) # Create parameter combination list 
		trials = 1000

		# main(grbp.name, template_grb, param_list, trials) # Run simulations

		plot(grbp.name, grbp.t_true) # Plot simulation results
	"""

	plot = PLOTSAMPLE()
	for i in range(len(grbs_names)):
		grbp = importlib.import_module("data_files.grb_{}.info".format(grbs_names[i]), package=None) # Load GRB parameters
		sim_results = np.load("data_files/grb_{}/grb_{}_redshift_sim-results.tmp.txt.npy".format(grbp.name, grbp.name))
		plot.add_data_table(sim_results)

	plot.cumulative_durations(keep_sep=False, z_min=3, label="Sim high z")
	ax = plt.gca()

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

	low_z_durations = np.zeros(shape=len(grbs_names),dtype=[("DURATION",float)])
	for i in range(len(grbs_names)):
		grbp = importlib.import_module("data_files.grb_{}.info".format(grbs_names[i]), package=None) # Load GRB parameters
		low_z_durations['DURATION'][i] = grbp.t_true

	plot.clear_data_table()
	plot.add_data_table(low_z_durations)
	plot.cumulative_durations(ax=ax, bin_max=np.log10(ax.get_xlim()[1]), label="Obs low z")

	high_z_durations = np.zeros(shape=len(high_z_grbs_names),dtype=[("DURATION",float)])
	for i in range(len(high_z_grbs_names)):
		grbp = importlib.import_module("data_files.grb_{}.info".format(high_z_grbs_names[i]), package=None) # Load GRB parameters
		high_z_durations['DURATION'][i] = grbp.t_true

	plot.clear_data_table()
	plot.add_data_table(high_z_durations)
	plot.cumulative_durations(ax=ax, bin_max=np.log10(ax.get_xlim()[1]), label="Obs high z")

	# plot.savefig(fname="data_files/figs/2024-03-27/cum_dur.png")
	plot.show()



