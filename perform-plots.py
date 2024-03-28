import importlib
import numpy as np 
import matplotlib.pyplot as plt
import time
from datetime import date

from packages.class_GRB import GRB
from packages.class_PLOTS import PLOTSIMRES, PLOTSAMPLE
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


def cum_dist_sep():
	plot = PLOTSAMPLE()
	ax = plt.figure().gca()

	for i in range(len(obs_low_z_grbs)):
		grbp = importlib.import_module("data_files.grb_{}.info".format(obs_low_z_grbs[i]), package=None) # Load GRB parameters
		sim_results = np.load("data_files/grb_{}/grb_{}_redshift_sim-results.tmp.txt.npy".format(grbp.name, grbp.name))
		sim_results = sim_results[sim_results['z'] <1]

		if len(sim_results) > 0:
			plot.cumulative_durations(data = sim_results, bin_max = 2e3, ax=ax)

	# plot.savefig(fname="data_files/figs/cum_dur_sep.png")

def cum_dist():
	plot = PLOTSAMPLE()

	sim_high_z = np.zeros(shape=0, dtype=dt_sim_res)
	for i in range(len(obs_low_z_grbs)):
		grbp = importlib.import_module("data_files.grb_{}.info".format(obs_low_z_grbs[i]), package=None) # Load GRB parameters
		sim_results = np.load("data_files/grb_{}/grb_{}_redshift_sim-results.tmp.txt.npy".format(grbp.name, grbp.name))
		sim_high_z = np.append(sim_high_z, sim_results[(sim_results["z"]>3) & (sim_results["z"]<9)])

	plot.cumulative_durations(data = sim_high_z, label="Sim high z")
	# plot.cumulative_fluence(data = sim_high_z, label="Sim high z")
	ax = plt.gca()

	sim_low_z = np.zeros(shape=0, dtype=dt_sim_res)
	for i in range(len(obs_low_z_grbs)):
		grbp = importlib.import_module("data_files.grb_{}.info".format(obs_low_z_grbs[i]), package=None) # Load GRB parameters
		sim_results = np.load("data_files/grb_{}/grb_{}_redshift_sim-results.tmp.txt.npy".format(grbp.name, grbp.name))
		sim_low_z = np.append(sim_low_z, sim_results[sim_results["z"]<1])
		# sim_low_z = np.append(sim_low_z, sim_results)

	plot.cumulative_durations(data = sim_low_z, ax=ax, bin_max=np.log10(ax.get_xlim()[1]), label="Sim low z")
	# plot.cumulative_fluence(data = sim_low_z, ax=ax, bin_max=np.log10(ax.get_xlim()[1]), label="Sim low z")

	obs_low_z = np.zeros(shape=len(obs_low_z_grbs),dtype=[("DURATION",float)])
	for i in range(len(obs_low_z_grbs)):
		grbp = importlib.import_module("data_files.grb_{}.info".format(obs_low_z_grbs[i]), package=None) # Load GRB parameters
		obs_low_z['DURATION'][i] = grbp.t_true

	plot.cumulative_durations(data=obs_low_z, ax=ax, bin_max=np.log10(ax.get_xlim()[1]), label="Obs low z")


	obs_high_z = np.zeros(shape=len(obs_high_z_grbs),dtype=[("DURATION",float)])
	for i in range(len(obs_high_z_grbs)):
		grbp = importlib.import_module("data_files.grb_{}.info".format(obs_high_z_grbs[i]), package=None) # Load GRB parameters
		obs_high_z['DURATION'][i] = grbp.t_true

	plot.cumulative_durations(data = obs_high_z, ax=ax, bin_max=np.log10(ax.get_xlim()[1]), label="Obs high z")

	# plot.savefig(fname="data_files/figs/cum_dur.png")

def cum_dist_fluence_sep():
	plot = PLOTSAMPLE()
	ax = plt.figure().gca()

	for i in range(len(obs_low_z_grbs)):
		grbp = importlib.import_module("data_files.grb_{}.info".format(obs_low_z_grbs[i]), package=None) # Load GRB parameters
		sim_results = np.load("data_files/grb_{}/grb_{}_redshift_sim-results.tmp.txt.npy".format(grbp.name, grbp.name))
		sim_results = sim_results[sim_results['z'] > 3]

		if len(sim_results) > 0:
			plot.cumulative_fluence(data = sim_results, bin_max = 2e3, ax=ax)

	# plot.savefig(fname="data_files/figs/cum_dur_sep.png")

def cum_dist_fluence():
	plot = PLOTSAMPLE()

	sim_high_z = np.zeros(shape=0, dtype=dt_sim_res)
	for i in range(len(obs_low_z_grbs)):
		grbp = importlib.import_module("data_files.grb_{}.info".format(obs_low_z_grbs[i]), package=None) # Load GRB parameters
		sim_results = np.load("data_files/grb_{}/grb_{}_redshift_sim-results.tmp.txt.npy".format(grbp.name, grbp.name))
		sim_high_z = np.append(sim_high_z, sim_results[(sim_results["z"]>3) & (sim_results["z"]<9)])

	plot.cumulative_fluence(data = sim_high_z, label="Sim high z")
	ax = plt.gca()

	sim_low_z = np.zeros(shape=0, dtype=dt_sim_res)
	for i in range(len(obs_low_z_grbs)):
		grbp = importlib.import_module("data_files.grb_{}.info".format(obs_low_z_grbs[i]), package=None) # Load GRB parameters
		sim_results = np.load("data_files/grb_{}/grb_{}_redshift_sim-results.tmp.txt.npy".format(grbp.name, grbp.name))
		sim_low_z = np.append(sim_low_z, sim_results[sim_results["z"]<1])
		# sim_low_z = np.append(sim_low_z, sim_results)

	plot.cumulative_fluence(data = sim_low_z, ax=ax, bin_max=np.log10(ax.get_xlim()[1]), label="Sim low z")

	# obs_low_z = np.zeros(shape=len(obs_low_z_grbs),dtype=[("DURATION",float)])
	# for i in range(len(obs_low_z_grbs)):
	# 	grbp = importlib.import_module("data_files.grb_{}.info".format(obs_low_z_grbs[i]), package=None) # Load GRB parameters
	# 	obs_low_z['DURATION'][i] = grbp.t_true

	# plot.cumulative_fluence(data=obs_low_z, ax=ax, bin_max=np.log10(ax.get_xlim()[1]), label="Obs low z")

	# obs_high_z = np.zeros(shape=len(obs_high_z_grbs),dtype=[("DURATION",float)])
	# for i in range(len(obs_high_z_grbs)):
	# 	grbp = importlib.import_module("data_files.grb_{}.info".format(obs_high_z_grbs[i]), package=None) # Load GRB parameters
	# 	obs_high_z['DURATION'][i] = grbp.t_true

	# plot.cumulative_fluence(data = obs_high_z, ax=ax, bin_max=np.log10(ax.get_xlim()[1]), label="Obs high z")

	# plot.savefig(fname="data_files/figs/cum_dur.png")

def z_evo():
	for i in range(len(obs_low_z_grbs)):
		grbp = importlib.import_module("data_files.grb_{}.info".format(obs_low_z_grbs[i]), package=None) # Load GRB parameters

		sim_results = np.load("data_files/grb_{}/grb_{}_redshift_sim-results.tmp.txt.npy".format(grbp.name, grbp.name))
		plot = PLOTSIMRES() # Plot simulation results
		plot.redshift_evo(sim_results, t_true=grbp.t_true, log=False)

		# plot.savefig(fname="data_files/figs/z-evo-plots/grb-{}-redshift-evo.png".format(grbp.name))


if __name__ == "__main__":

	# cum_dist_sep()
	# cum_dist()
	cum_dist_fluence_sep()
	# cum_dist_fluence()
	# z_evo() # z_evo plots for each GRB

	plt.show()
