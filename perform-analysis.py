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


# Make a GRB object
z = 1  # redshift 
template_grb = GRB(z=z)
# Make light curve 
template_grb.load_light_curve("data-files/template-light-curves/grb_110422A_1chan_1s.lc", rm_trigtime=True, det_area=0.16)
template_grb.light_curve = template_grb.light_curve[np.argmax(-100 <= template_grb.light_curve['TIME']):np.argmax(template_grb.light_curve['TIME'] >= 200)]

alpha = -0.831 
ep = 147.995
tmin = -10
tmax = 35
norm = template_grb.get_ave_photon_flux(tmin=tmin,tmax=tmax) / CPL(alpha= alpha,ep=ep,norm=1,enorm=1)._calc_phot_flux(15, 150) # for GRB GRB 081007
template_grb.load_specfunc(CPL(alpha= alpha,ep=ep,norm=norm,enorm=1))


# Simulate many observations 
z_arr = np.array([1])

imx_arr = np.array([0.,0.8,1.2])
# imx_arr = np.array([0.])
# imx_arr = np.linspace(-1,1,num=10)

imy_arr = np.array([0.])
# imy_arr = np.linspace(-0.7,0.7,num=10)

ndets_arr = np.array([30000])

param_list = make_param_list(z_arr,imx_arr,imy_arr,ndets_arr)
trials = 3




sim_results, grbs = many_simulations(template_grb, param_list, trials, multiproc=False, ret_ave=True, keep_synth_grbs=True)


template_grb.light_curve['RATE'] *= 0.16



plots = PLOTS(sim_results)
# plots.det_plane_map()
# plots.dur_vs_param(obs_param="z",t_true=12,dur_frac=True)
plots.plot_light_curves(grbs=np.append(template_grb,grbs),labels=["Template",*imx_arr],t_window=[-40,60])
plt.show()
