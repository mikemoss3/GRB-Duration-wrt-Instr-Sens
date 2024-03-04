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
name = "060614"
z = 0.1254
fn = "data-files/grb-060614/grb_060614_1chan_64ms.lc"
t_true = 109.104 # true T90
alpha = -2.33
ep = 134.112
norm = 9.47830E-03 # counts cm−2 s^−1 keV^−1


# Make a GRB object
template_grb = GRB(z=z)
# Load light curve 
template_grb.load_light_curve(fn, rm_trigtime=True, det_area=0.16)
template_grb.cut_light_curve(tmin=-70, tmax=230)

# Load spectrum
template_grb.load_specfunc(CPL(alpha= alpha, ep=ep, norm=norm, enorm=50))

# z_arr = np.array([z])
z_arr = np.linspace(z, 1.3, num=50)

imx_arr = np.array([0.])
imy_arr = np.array([0.])
# imx_arr = np.linspace(-1.75,1.75,70)
# imy_arr = np.linspace(-0.875,0.875,70)

ndets_arr = np.array([30000])

param_list = make_param_list(z_arr,imx_arr,imy_arr,ndets_arr)
trials = 100

def main():
	sim_results = many_simulations(template_grb, param_list, trials, multiproc=False, keep_synth_grbs=False, bgd_rate_per_det=0.3, verbose=True)
	ave_sim_results = make_ave_sim_res(sim_results)

	np.save("data-files/grb-{}/grb_{}_redshift_sim-results.tmp.txt".format(name, name), sim_results)
	np.save("data-files/grb-{}/grb_{}_redshift_ave-sim-results.tmp.txt".format(name, name), ave_sim_results)

	# np.save("data-files/grb-{}/grb_{}_detector_sim-results.tmp.txt".format(name, name), sim_results)
	# np.save("data-files/grb-{}/grb_{}_detector_ave-sim-results.tmp.txt".format(name, name), ave_sim_results)
	
	# for i in range(len(grbs)):
		# np.savetxt("/Users/mjmoss/Research/presentation-plot-making/2024-02-15-SED-seminar/grb-{}/light_curve_at_{}.txt".format(name, z_arr[i]), grbs[i].light_curve)


def plot():
	sim_results = np.load("data-files/grb-{}/grb_{}_redshift_sim-results.tmp.txt.npy".format(name, name))
	plot = PLOTS()
	plot.redshift_evo(sim_results, t_true=t_true, log=False)
	

	# ave_sim_results = np.load("data-files/grb-{}/grb_{}_detector_ave-sim-results.tmp.txt.npy".format(name, name))
	# plot = PLOTS()
	# plot.det_plane_map(ave_sim_results, inc_grids=True)


	# today = date.today()
	# plt.savefig("data-files/figs/{}/{}-grb-{}-redshift-evo.png".format(today, today, name), dpi=400)
	
	plt.show()


if __name__ == "__main__":
	# main()
	plot()
