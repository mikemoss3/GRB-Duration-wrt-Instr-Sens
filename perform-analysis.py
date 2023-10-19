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

z = 1  # redshift 
imx, imy = 0., 0.  # Position on the detector plane
ndets = 30000  # Number of enabled detectors

# Make a GRB object
template_grb = GRB(z=z)
# Make light curve 
template_grb.load_light_curve("data-files/template-light-curves/grb_160121A_1chan_1s.lc", rm_trigtime=True)
template_grb.light_curve = template_grb.light_curve[np.argmax(-100 <= template_grb.light_curve['TIME']):np.argmax(template_grb.light_curve['TIME'] >= 100)]

template_grb.load_specfunc(PL(alpha=-1.77,norm=4.43e-3,enorm=50)) # for GRB GRB 081007

# Simulate many observations 
z_arr = np.array([1])
imx_arr = np.array([0.,0.5,0.8])
# imx_arr = np.linspace(-1,1,num=10)
imy_arr = np.array([0.])
# imy_arr = np.linspace(-0.7,0.7,num=10)
ndets_arr = np.array([30000])
param_list = make_param_list(z_arr,imx_arr,imy_arr,ndets_arr)
trials = 10

sim_results, grbs = many_simulations(template_grb, param_list, trials, multiproc=False, ret_ave=True, keep_synth_grbs=True)

plots = PLOTS(sim_results)
# plots.det_plane_map()
# plots.dur_vs_param(obs_param="z",t_true=12,dur_frac=True)
plots.plot_light_curves(grbs=np.append(template_grb,grbs),labels=["Template",*imx_arr])
plt.show()

