"""
Author: Mike Moss
Contact: mikejmoss3@gmail.com	

Sandbox to perform all simulation, analysis, and plotting commands

"""


import numpy as np 
import matplotlib.pyplot as plt
from datetime import date

from packages.class_GRB import GRB
from packages.class_PLOTS import PLOTS
from packages.class_SPECFUNC import PL, CPL
from packages.package_many_simulations import many_simulations, make_param_list, make_ave_sim_res
from util_packages.package_datatypes import dt_sim_res


# GRB info
name = "080810"
z = 3.3604
fn = "data-files/grb-080810/grb_080810_1chan_1s.lc"
t_true = 107.668 # true T90
alpha = -1.270
ep = 9999.350
norm = 1.39e-02 # counts cm−2 s^−1 keV^−1



# Make a GRB object
template_grb = GRB(z=z)
# Load light curve 
template_grb.load_light_curve(fn, rm_trigtime=True, det_area=0.16)
template_grb.cut_light_curve(tmin=-100, tmax=250)
# template_grb.zero_light_curve_selection(tmax=tmin)
# template_grb.zero_light_curve_selection(tmin=tmax)

# Load spectrum
template_grb.load_specfunc(CPL(alpha= alpha, ep=ep, norm=norm, enorm=50))

z_arr = np.array([z])
z_arr = np.linspace(z, 12, num=50)

imx_arr = np.array([0])
imy_arr = np.array([0])
# imx_arr = np.linspace(-1.75,1.75,75)
# imy_arr = np.linspace(-0.875,0.875,75)

ndets_arr = np.array([30000])

param_list = make_param_list(z_arr,imx_arr,imy_arr,ndets_arr)
trials = 1000

def main():
	sim_results = many_simulations(template_grb, param_list, trials, multiproc=False, keep_synth_grbs=False, bgd_rate_per_det=0.3, verbose=True)
	ave_sim_results = make_ave_sim_res(sim_results)

	np.save("data-files/grb-{}/grb_{}_redshift_sim-results.txt".format(name, name), sim_results)
	np.save("data-files/grb-{}/grb_{}_redshift_ave-sim-results.txt".format(name, name), ave_sim_results)

	# np.save("data-files/grb-{}/grb_{}_detector_sim-results.txt".format(name, name), sim_results)
	# np.save("data-files/grb-{}/grb_{}_detector_ave-sim-results.txt".format(name, name), ave_sim_results)
	
	# template_grb.light_curve['RATE'] *= 0.16
	# plt = PLOTS()
	# plt.plot_light_curves(grbs=np.append(template_grb,grbs),t_window=[-50, 150])
	# plt.show()


def plot():

	sim_results = np.load("data-files/grb-{}/grb_{}_redshift_sim-results.txt.npy".format(name, name))
	ave_sim_results = np.load("data-files/grb-{}/grb_{}_redshift_ave-sim-results.txt.npy".format(name, name))
	plt = PLOTS()
	plt.redshift_evo(sim_results, t_true=t_true)
	
	
	"""
	sim_results = np.load("data-files/grb-{}/grb_{}_detector_sim-results.txt.npy".format(name, name))
	ave_sim_results = np.load("data-files/grb-{}/grb_{}_detector_ave-sim-results.txt.npy".format(name, name))
	plt = PLOTS()
	plt.det_plane_map(ave_sim_results, inc_grids=True)
	"""

	# today = date.today()
	# plt.savefig("data-files/figs/{}/{}-grb-{}.png".format(today, today, name))
	
	plt.show()


if __name__ == "__main__":
	# main()
	plot()
