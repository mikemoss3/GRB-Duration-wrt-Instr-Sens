import numpy as np 
import matplotlib.pyplot as plt

from packages.class_GRB import GRB
from packages.class_RSP import ResponseMatrix
from packages.class_SPECFUNC import PL, CPL
from packages.package_simulations import simulate_observation
from packages.package_analysis import many_simulations


z = 1  # redshift 
imx, imy = 0., 0.  # Position on the detector plane
ndets = 30000  # Number of enabled detectors

# Make a GRB object
template_grb = GRB(z=z,imx=imx,imy=imy)
# Make light curve 
template_grb.load_light_curve("data-files/sw00634795000b_1chan_1s.lc", rm_trigtime=True)
template_grb.light_curve = template_grb.light_curve[np.argmax(-200 <= template_grb.light_curve['TIME']):np.argmax(template_grb.light_curve['TIME'] >= 200)]

template_grb.load_specfunc(CPL(alpha=-1.,norm=1))

# Make a Response Matrix object and Load Swift BAT response based on position on detector plane  
resp_mat = ResponseMatrix()
resp_mat.load_SwiftBAT_resp(imx, imy)

# Simulate many observations 
z_arr = np.array([1])
imx_arr, imy_arr = np.array([0]), np.array([0])
ndets_arr = np.array([30000])
trials = 10
sim_results = many_simulations(template_grb, resp_mat, z_arr, imx_arr, imy_arr, ndets_arr, trials,multiproc=False)

print(sim_results)

# Plot results 

