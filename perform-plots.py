import importlib
import numpy as np 
import matplotlib.pyplot as plt
import time
from datetime import date

from packages.class_GRB import GRB
from packages.class_PLOTS import PLOTSIMRES, PLOTSAMPLE
from util_packages.package_datatypes import dt_sim_res
from packages.package_many_simulations import make_ave_sim_res


obs_low_z_grbs = np.array([
	# "050416A",
	# "050525A",
	# "060614",
	# "060912A",
	# "061021",
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
	"161219B",
	], dtype="U11")

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
		sim_results = np.load("data_files/results_final/grb_{}_redshift_sim-results.txt.npy".format(grbp.name, grbp.name))

		if len(sim_results) > 0:
			# plot.cumulative_durations(data = sim_results, bin_max = 2e3, ax=ax)
			plot.cumulative_durations(data = sim_results, bin_min=0.01, bin_max = 1, ax=ax, normed=True)

	# plot.savefig(fname="data_files/figs/2024-04-09/cum_dur_sep_max-normed.png")

def cum_dist():
	plot = PLOTSAMPLE()

	sim_high_z = np.zeros(shape=0, dtype=dt_sim_res)
	for i in range(len(obs_low_z_grbs)):
		grbp = importlib.import_module("data_files.grb_{}.info".format(obs_low_z_grbs[i]), package=None) # Load GRB parameters
		sim_results = np.load("data_files/results_final/grb_{}_redshift_sim-results.txt.npy".format(grbp.name, grbp.name))
		sim_high_z = np.append(sim_high_z, sim_results[(sim_results["z"]>3) & (sim_results["z"]<9)])

	plot.cumulative_durations(data = sim_high_z, label="Sim high z")
	ax = plt.gca()

	sim_low_z = np.zeros(shape=0, dtype=dt_sim_res)
	for i in range(len(obs_low_z_grbs)):
		grbp = importlib.import_module("data_files.grb_{}.info".format(obs_low_z_grbs[i]), package=None) # Load GRB parameters
		sim_results = np.load("data_files/results_final/grb_{}_redshift_sim-results.txt.npy".format(grbp.name, grbp.name))
		sim_low_z = np.append(sim_low_z, sim_results[sim_results["z"]<1])

	# plot.cumulative_durations(data = sim_low_z, ax=ax, bin_max=ax.get_xlim()[1], label="Sim low z")

	obs_low_z = np.zeros(shape=len(obs_low_z_grbs),dtype=[("DURATION",float)])
	for i in range(len(obs_low_z_grbs)):
		grbp = importlib.import_module("data_files.grb_{}.info".format(obs_low_z_grbs[i]), package=None) # Load GRB parameters
		obs_low_z['DURATION'][i] = grbp.t_true

	plot.cumulative_durations(data=obs_low_z, ax=ax, bin_max=ax.get_xlim()[1], label="Obs low z")

	obs_high_z = np.zeros(shape=len(obs_high_z_grbs),dtype=[("z",float),("DURATION",float)])
	for i in range(len(obs_high_z_grbs)):
		grbp = importlib.import_module("data_files.grb_{}.info".format(obs_high_z_grbs[i]), package=None) # Load GRB parameters
		obs_high_z['z'][i] = grbp.z
		obs_high_z['DURATION'][i] = grbp.t_true

	plot.cumulative_durations(data = obs_high_z, ax=ax, bin_max=ax.get_xlim()[1], label="Obs high z")

	# plot.savefig(fname="data_files/figs/2024-04-09/cum_dur.png")

def cum_dist_fluence_sep():
	plot = PLOTSAMPLE()
	ax = plt.figure().gca()

	for i in range(len(obs_low_z_grbs)):
		grbp = importlib.import_module("data_files.grb_{}.info".format(obs_low_z_grbs[i]), package=None) # Load GRB parameters
		sim_results = np.load("data_files/results_final/grb_{}_redshift_sim-results.txt.npy".format(grbp.name, grbp.name))

		if len(sim_results) > 0:
			# sim_results['FLUENCE']/=np.min(sim_results['FLUENCE'][sim_results['FLUENCE']>0])
			plot.cumulative_fluence(data = sim_results, bin_max = 3e3, ax=ax)

	# plot.savefig(fname="data_files/figs/2024-04-09/cum_fluence_sep.png")

def cum_dist_peak_flux_sep():
	plot = PLOTSAMPLE()
	ax = plt.figure().gca()

	for i in range(len(obs_low_z_grbs)):
		grbp = importlib.import_module("data_files.grb_{}.info".format(obs_low_z_grbs[i]), package=None) # Load GRB parameters
		sim_results = np.load("data_files/results_final/grb_{}_redshift_sim-results.txt.npy".format(grbp.name, grbp.name))

		if len(sim_results) > 0:
			# sim_results['1sPeakFlux']/=np.min(sim_results['1sPeakFlux'][sim_results['1sPeakFlux']>0])
			plot.cumulative_peak_flux(data = sim_results, bin_max = 1e3, ax=ax)

	ax.axvline(x=1.18 * 2.4*10**(-2) * 0.064**-0.5,ymin=0,ymax=1, linestyle="dashed", color="k")
	ax.axvline(x=1.18 * 2.4*10**(-2) * 1**-0.5,ymin=0,ymax=1, linestyle="dashed", color="purple")

	# plot.savefig(fname="data_files/figs/2024-04-09/cum_peak_flux_sep_v01.png")

def cum_dist_peak_flux():
	plot = PLOTSAMPLE()

	sim_high_z = np.zeros(shape=0, dtype=dt_sim_res)
	sim_high_z_1s = np.zeros(shape=0, dtype=dt_sim_res)
	sim_high_z_64ms = np.zeros(shape=0, dtype=dt_sim_res)
	for i in range(len(obs_low_z_grbs)):
		grbp = importlib.import_module("data_files.grb_{}.info".format(obs_low_z_grbs[i]), package=None) # Load GRB parameters
		sim_results = np.load("data_files/results_final/grb_{}_redshift_sim-results.txt.npy".format(grbp.name, grbp.name))

		sim_high_z = np.append(sim_high_z, sim_results[(sim_results["z"]>3) & (sim_results["z"]<9)])
		
		# if grbp.fn[-5:] == "1s.lc":
		# 	sim_high_z_1s = np.append(sim_high_z_1s, sim_results[(sim_results["z"]>3) & (sim_results["z"]<9)])
		# elif grbp.fn[-7:] == "64ms.lc":
		# 	sim_high_z_64ms = np.append(sim_high_z_64ms, sim_results[(sim_results["z"]>3) & (sim_results["z"]<9)])

	# sim_high_z_64ms['1sPeakFlux']*=0.064 
	plot.cumulative_peak_flux(data = sim_high_z, bin_max=10, label="Sim high z")
	ax = plt.gca()
	# plot.cumulative_peak_flux(data = sim_high_z_1s, ax=ax, bin_max=10, label="Sim high z - 1s")
	# plot.cumulative_peak_flux(data = sim_high_z_64ms, ax=ax, bin_max=10, label="Sim high z - 64ms")

	sim_low_z = np.zeros(shape=0, dtype=dt_sim_res)
	sim_low_z_1s = np.zeros(shape=0, dtype=dt_sim_res)
	sim_low_z_64ms = np.zeros(shape=0, dtype=dt_sim_res)
	for i in range(len(obs_low_z_grbs)):
		grbp = importlib.import_module("data_files.grb_{}.info".format(obs_low_z_grbs[i]), package=None) # Load GRB parameters
		sim_results = np.load("data_files/results_final/grb_{}_redshift_sim-results.txt.npy".format(grbp.name, grbp.name))
		sim_low_z = np.append(sim_low_z, sim_results[sim_results["z"]<1])

		# if grbp.fn[-5:] == "1s.lc":
		# 	sim_low_z_1s = np.append(sim_low_z_1s, sim_results[sim_results["z"]<1])
		# elif grbp.fn[-7:] == "64ms.lc":
		# 	sim_low_z_64ms = np.append(sim_low_z_64ms, sim_results[sim_results["z"]<1])

	plot.cumulative_peak_flux(data = sim_low_z, ax=ax, bin_max=ax.get_xlim()[1], label="Sim low z")
	# plot.cumulative_peak_flux(data = sim_low_z_1s, ax=ax, bin_max=ax.get_xlim()[1], label="Sim low z - 1s")
	# plot.cumulative_peak_flux(data = sim_low_z_64ms, ax=ax, bin_max=ax.get_xlim()[1], label="Sim low z - 64ms")

	# obs_low_z = np.zeros(shape=len(obs_low_z_grbs),dtype=[("DURATION",float)])
	# for i in range(len(obs_low_z_grbs)):
	# 	grbp = importlib.import_module("data_files.grb_{}.info".format(obs_low_z_grbs[i]), package=None) # Load GRB parameters
	# 	obs_low_z['DURATION'][i] = grbp.t_true

	# plot.cumulative_peak_flux(data=obs_low_z, ax=ax, bin_max=ax.get_xlim()[1], label="Obs low z")

	obs_high_z_peak_fluxes = np.array([0.3696, 0.167, 0.7110, 0.356, 2.927, 0.487, 0.808, 0.298, 0.740, 0.3409, 0.6056,])

	obs_high_z = np.zeros(shape=len(obs_high_z_grbs),dtype=[("1sPeakFlux",float),("z",float)])
	for i in range(len(obs_high_z_grbs)):
		grbp = importlib.import_module("data_files.grb_{}.info".format(obs_high_z_grbs[i]), package=None) # Load GRB parameters
		obs_high_z['z'][i] = grbp.z
		obs_high_z['1sPeakFlux'][i] = obs_high_z_peak_fluxes[i]

	plot.cumulative_peak_flux(data = obs_high_z, ax=ax, bin_max=10, label="Obs high z")

	ax.axvline(x=1.18 * 2.4*10**(-2) * 0.064**-0.5,ymin=0,ymax=1, linestyle="dashed", color="k")
	ax.axvline(x=1.18 * 2.4*10**(-2) * 1**-0.5,ymin=0,ymax=1, linestyle="dashed", color="purple")

	# plot.savefig(fname="data_files/figs/2024-04-09/cum_peak_flux_v01.png")

def redshift_dist():
	plot = PLOTSAMPLE()

	sim_high_z = np.zeros(shape=0, dtype=dt_sim_res)
	z_array = np.zeros(shape=12, dtype=[("z_min",float),("z_max",float),("numGRBs",int)])
	z_array['z_min'] = np.linspace(3, 10, num=len(z_array))
	z_array['z_max'] = np.linspace(3.5, 10.5, num=len(z_array))
	for i in range(len(obs_low_z_grbs)):
		grbp = importlib.import_module("data_files.grb_{}.info".format(obs_low_z_grbs[i]), package=None) # Load GRB parameters
		sim_results = np.load("data_files/results_final/grb_{}_redshift_sim-results.txt.npy".format(grbp.name, grbp.name))
		unique_z = np.unique(sim_results['z'])
		for j in range(len(z_array)):
			if any(unique_z > z_array['z_min'][j]):
				z_array['numGRBs'][j] += 1

	ax = plt.figure().gca()
	ax.step(x=z_array['z_min'], y=z_array['numGRBs'], where="post")
	plot.plot_aesthetics(ax)
	ax.margins(x=0.1,y=0.1)
	ax.set_ylim(0)
	ax.set_xlim(3)
	from matplotlib.ticker import MaxNLocator
	ax.yaxis.set_major_locator(MaxNLocator(integer=True))
	ax.set_xlabel("z", fontsize=14)
	ax.set_ylabel("Unique GRBs Still Detected", fontsize=14)

def z_evo():
	for i in range(len(obs_low_z_grbs)):
		grbp = importlib.import_module("data_files.grb_{}.info".format(obs_low_z_grbs[i]), package=None) # Load GRB parameters

		# sim_results = np.load("data_files/results_final/grb_{}_redshift_sim-results.txt.npy".format(grbp.name, grbp.name))
		sim_results = np.load("data_files/grb_{}/grb_{}_redshift_sim-results.tmp.txt.npy".format(grbp.name, grbp.name, grbp.name))
		plot = PLOTSIMRES() # Plot simulation results
		plot.redshift_evo(sim_results, t_true=grbp.t_true, log=False)

		# plot.savefig(fname="data_files/figs/z-evo-plots/grb-{}-redshift-evo.png".format(grbp.name), dpi="figure")
		# plot.close()

def z_fluence_evo():
	for i in range(len(obs_low_z_grbs)):
		grbp = importlib.import_module("data_files.grb_{}.info".format(obs_low_z_grbs[i]), package=None) # Load GRB parameters

		# sim_results = np.load("data_files/results_final/grb_{}_redshift_sim-results.txt.npy".format(grbp.name, grbp.name))
		sim_results = np.load("data_files/grb_{}/grb_{}_redshift_sim-results.tmp.txt.npy".format(grbp.name, grbp.name, grbp.name))
		plot = PLOTSIMRES() # Plot simulation results
		plot.redshift_fluence_evo(sim_results)
		
		# plot.savefig(fname="data_files/figs/z-evo-fluence-plots/grb-{}-redshift-evo.png".format(grbp.name), dpi="figure")
		# plot.close()

def t90_vs_z_below_2s():
	
	plot = PLOTSIMRES()
	ax = plt.gca()
	
	for j in range(len(obs_low_z_grbs)):
		grbp = importlib.import_module("data_files.grb_{}.info".format(obs_low_z_grbs[j]), package=None) # Load GRB parameters
		sim_results = np.load("data_files/results_final/grb_{}_redshift_sim-results.txt.npy".format(grbp.name, grbp.name))

		ave_sim_results = make_ave_sim_res(sim_results, omit_nondetections=True)

		if any(ave_sim_results['DURATION']<3):
			line = plot.dur_vs_param(ave_sim_results, 'z', ax=ax, label=grbp.name, linewidth=1.5, joined=True)

		"""
		# To use minimum instead of average T90
		uniq_z = np.unique(sim_results['z'])
		min_t90 = np.zeros(shape=len(uniq_z), dtype=[("z",float),("DURATION",float)])
		min_t90['z'] = uniq_z
		for i in range(len(uniq_z)):
			min_t90['DURATION'][i] = np.min(sim_results['DURATION'][sim_results['z']==uniq_z[i]])

		if any(min_t90['DURATION']<3):
			plot.dur_vs_param(min_t90, 'z', ax=ax, label=grbp.name, linestyle="dashed", alpha=0.6, joined=True, color=line.get_color())
		"""

	ax.axhline(y=2,xmin=0,xmax=1, color="k", linestyle="dashed")

	# ax.set_yscale("log")
	# ax.set_ylim(0.06)
	ax.set_ylim(0, 20)
	ax.legend()

	ax.set_ylabel(r"Measured T$_{90}$", fontsize=14)
	ax.set_xlabel("z", fontsize=14)

	plot.plot_aesthetics(ax)


if __name__ == "__main__":

	# cum_dist_sep()
	# cum_dist()
	# cum_dist_fluence_sep()
	# cum_dist_peak_flux_sep()
	# cum_dist_peak_flux()
	# redshift_dist()
	z_evo()
	z_fluence_evo()
	# t90_vs_z_below_2s()

	plt.show()
