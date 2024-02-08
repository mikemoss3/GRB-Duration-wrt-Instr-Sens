"""
Author: Mike Moss
Contact: mikejmoss3@gmail.com	

Sandbox to perform all simulation, analysis, and plotting commands

"""


import numpy as np 
import matplotlib.pyplot as plt

from packages.class_GRB import GRB
from packages.class_PLOTS import PLOTS
from packages.class_SPECFUNC import PL, CPL
from packages.package_many_simulations import many_simulations, make_param_list, make_ave_sim_res


def main():
	# GRB info
	z = 1.758
	fn = "data-files/template-light-curves/grb_150314A_1chan_1s.lc"
	t_true = 14.780 # true T90
	tmin = -20
	tmax = 170
	alpha = -0.465417
	ep = 100.872
	norm = 0.01 # counts cm−2 s^−1 keV^−1

	# Make a GRB object
	template_grb = GRB(z=z)
	# Load light curve 
	template_grb.load_light_curve(fn, rm_trigtime=True, det_area=0.16)
	template_grb.light_curve = template_grb.light_curve[np.argmax(tmin <= template_grb.light_curve['TIME']):np.argmax(template_grb.light_curve['TIME'] >= tmax)]

	# Load spectrum
	template_grb.load_specfunc(CPL(alpha= alpha, ep=ep, norm=norm, enorm=50))



	# Simulate many observations 
	z_arr = np.array([z])
	# z_arr = np.linspace(z, 4, num=4)

	imx_arr = np.array([0])
	imy_arr = np.array([0])

	imx_arr = np.linspace(-1.75,1.75,10)
	imy_arr = np.linspace(-0.875,0.875,10)

	ndets_arr = np.array([20000])

	param_list = make_param_list(z_arr,imx_arr,imy_arr,ndets_arr)
	trials = 4

	sim_results, grbs = many_simulations(template_grb, param_list, trials, multiproc=False, ret_ave=False, keep_synth_grbs=True, bgd_rate_per_det=0.5)
	ave_sim_results = make_ave_sim_res(sim_results)

	# print(sim_results)
	# print(ave_sim_results)

	template_grb.light_curve['RATE'] *= 0.16


	# plt = PLOTS()
	# plt.plot_light_curves(template_grb)
	# plt.plot_light_curves(grbs=np.append(template_grb,grbs),t_window=[tmin,tmax])
	# plt.plot_spectra(np.append(template_grb,grbs), emin=15, emax = 150, labels= np.append(["template"],z_arr))
	# plt.dur_vs_param(obs_param="z")
	# plt.redshift_evo(t_true=t_true)
	plt = PLOTS()
	plt.det_plane_map(ave_sim_results, inc_grids=True)
	
	plt.show()

	# plt.savefig("data-files/figs/2023-11-06/2023-11-06-grb-211211A-redshift-vs-dur-evo.png")

if __name__ == "__main__":
	main()
