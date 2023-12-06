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
	z = 0.5295
	fn = "data-files/template-light-curves/grb_081007_1chan_64ms.lc"
	t_true = 9.728 # true T90
	tmin = -15
	tmax = 20
	alpha = -1.52043
	ep = 20.1229
	norm = 1.77485E-02

	# Make a GRB object
	template_grb = GRB(z=z)
	# Load light curve 
	template_grb.load_light_curve(fn, rm_trigtime=True, det_area=0.16)
	template_grb.light_curve = template_grb.light_curve[np.argmax(tmin <= template_grb.light_curve['TIME']):np.argmax(template_grb.light_curve['TIME'] >= tmax)]

	# Load spectrum
	template_grb.load_specfunc(CPL(alpha= alpha, ep=ep, norm=norm, enorm=50))

	# Simulate many observations 
	# z_arr = np.array([0.1254,0.2,0.3]) # 060614
	# z_arr = np.array([0.076,0.3,0.5,0.7,0.9,1.1,1.3,1.4,1.5,1.6,1.7,2]) # 211211A

	# z_arr = np.linspace(z,8,50)
	z_arr = np.array([z])

	# imx_arr = np.array([0.])
	# imy_arr = np.array([0.])

	imx_arr = np.linspace(-1,1.,30)
	imy_arr = np.linspace(-0.7,0.7,30)

	ndets_arr = np.array([30000])

	param_list = make_param_list(z_arr,imx_arr,imy_arr,ndets_arr)
	trials = 1


	sim_results, grbs = many_simulations(template_grb, param_list, trials, multiproc=False, ret_ave=True, keep_synth_grbs=True,bgd_rate_per_det=0.6)


	template_grb.light_curve['RATE'] *= 0.16



	plots = PLOTS(sim_results)
	# plots.plot_light_curves(template_grb)
	plots.det_plane_map()
	# plots.dur_vs_param(obs_param="z",t_true=t_true,dur_frac=False)
	# plots.plot_light_curves(grbs=np.append(template_grb,grbs),labels=["Template",*z_arr],t_window=[tmin,tmax])
	plt.show()

	# plt.savefig("data-files/figs/2023-11-06/2023-11-06-grb-211211A-redshift-vs-dur-evo.png")

if __name__ == "__main__":
	main()
