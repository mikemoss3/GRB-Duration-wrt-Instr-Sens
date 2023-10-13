"""
Author: Mike Moss
Contact: mikejmoss3@gmail.com	

Sandbox to perform all simulation, analysis, and plotting commands

"""


import numpy as np 
import matplotlib.pyplot as plt

from packages.class_GRB import GRB
from packages.class_RSP import ResponseMatrix
from packages.class_SPECFUNC import PL, CPL
from packages.package_simulations import simulate_observation
from packages.package_analysis import many_simulations, make_ave_sim_res
from packages.package_plotting import PLOTS

z = 1  # redshift 
imx, imy = 0., 0.  # Position on the detector plane
ndets = 30000  # Number of enabled detectors

# Make a GRB object
template_grb = GRB(z=z,imx=imx,imy=imy)
# Make light curve 
# template_grb.load_light_curve("data-files/template-light-curves/grb_150314A_1chan_1s.lc", rm_trigtime=True)
template_grb.load_light_curve("data-files/template-light-curves/grb_180404A_1chan_1s.lc", rm_trigtime=True)
template_grb.light_curve = template_grb.light_curve[np.argmax(-200 <= template_grb.light_curve['TIME']):np.argmax(template_grb.light_curve['TIME'] >= 200)]

# template_grb.load_specfunc(CPL(alpha=-1.,norm=4)) # for GRB 150314A
template_grb.load_specfunc(CPL(alpha=-1.9,norm=60)) # for GRB 180404A

# Make a Response Matrix object and Load Swift BAT response based on position on detector plane  
resp_mat = ResponseMatrix()
resp_mat.load_SwiftBAT_resp(imx, imy)

# Simulate many observations 
z_arr = np.array([1])
imx_arr, imy_arr = np.array([0.]), np.array([0.])
ndets_arr = np.array([30000])
trials = 10
sim_results, grbs = many_simulations(template_grb, resp_mat, z_arr, imx_arr, imy_arr, ndets_arr, trials,multiproc=False,ret_ave=False,keep_synth_grbs=True)

print(sim_results)

plots = PLOTS(sim_results)
plots.plot_light_curves(grbs=np.append(template_grb,grbs),labels=["Template","z=1","z=2"])
plt.show()


