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
from packages.package_analysis import many_simulations, make_param_list


def main():
	# GRB info
	z = 1.758
	fn = "data-files/template-light-curves/grb_150314A_1chan_1s.lc"
	t_true = 14.780 # true T90
	tmin = -20
	tmax = 170
	alpha = -0.465417
	ep = 100.872
	norm = 0.03 # counts cm−2 s^−1 keV^−1

	# Make a GRB object
	template_grb = GRB(z=z)
	# Load light curve 
	template_grb.load_light_curve(fn, rm_trigtime=True, det_area=0.16)
	template_grb.light_curve = template_grb.light_curve[np.argmax(tmin <= template_grb.light_curve['TIME']):np.argmax(template_grb.light_curve['TIME'] >= tmax)]

	# Load spectrum
	template_grb.load_specfunc(CPL(alpha= alpha, ep=ep, norm=norm, enorm=50))



	# Simulate many observations 
	# z_arr = np.array([z])
	z_arr = np.linspace(z, 8, num=3)

	imx_arr = np.array([0])
	imy_arr = np.array([0])

	# imx_arr = np.linspace(-1.25,1.25,30)
	# imy_arr = np.linspace(-0.875,0.875,30)

	ndets_arr = np.array([30000])

	param_list = make_param_list(z_arr,imx_arr,imy_arr,ndets_arr)
	trials = 4

	sim_results = many_simulations(template_grb, param_list, trials, multiproc=False, ret_ave=False, keep_synth_grbs=False, bgd_rate_per_det=0.6)


	# print(sim_results)

	template_grb.light_curve['RATE'] *= 0.16



	plots = PLOTS(sim_results)
	# plots.plot_light_curves(template_grb)
	# plots.plot_light_curves(grbs=np.append(template_grb,grbs),t_window=[tmin,tmax])
	# plots.dur_vs_param(obs_param="z")
	# plots.det_plane_map()
	plots.redshift_evo(t_true=t_true)
	plt.show()

	# plt.savefig("data-files/figs/2023-11-06/2023-11-06-grb-211211A-redshift-vs-dur-evo.png")

if __name__ == "__main__":
	main()
