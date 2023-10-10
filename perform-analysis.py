import numpy as np 
import matplotlib.pyplot as plt

from packages.class_GRB import GRB
from packages.class_RSP import ResponseMatrix
from packages.package_simulations import simulate_observation
from packages.package_bayesian_block import bayesian_t_blocks

z = 1 # redshift 
imx, imy = 0., 0. # Position on the detector plane
ndets = 10000 # Number of enabled detectors

# Make a GRB object
template_grb = GRB()
# Make light curve 
template_grb.load_light_curve("data-files/sw00330856000b_1chan_64ms.lc",rm_trigtime=True)
template_grb.light_curve = template_grb.light_curve[np.argmax(-200<=template_grb.light_curve['TIME']):np.argmax(template_grb.light_curve['TIME']>=200)]

template_grb.load_spectrum("data-files/test_spectrum.txt")

template_grb.move_to_source_frame(z=z)

# Make a Response Matrix object 
resp_mat = ResponseMatrix()
# Load Swift BAT response based on position on detector plane 
resp_mat.load_SwiftBAT_resp(imx, imy)

# Simulate observation
synth_grb = simulate_observation(template_grb, z, imx, imy, ndets, resp_mat)


# Measure T90 using Bayesian block
synth_t90, t_start, t_end, synth_t90_flue = bayesian_t_blocks(synth_grb,dur_per=90)

print(synth_t90, t_start, t_end, synth_t90_flue)

plt.plot(template_grb.light_curve['TIME'],template_grb.light_curve['RATE'])
plt.plot(synth_grb.light_curve['TIME'],synth_grb.light_curve['RATE'])
# plt.axvspan(xmin=-1.2,xmax=13.5,color="C0",alpha=0.3)
plt.axvspan(xmin=t_start,xmax=t_end,color="C3",alpha=0.3)

plt.xlim(-75,75)

plt.show()

