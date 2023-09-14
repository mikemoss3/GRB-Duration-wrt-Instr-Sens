import numpy as np 

from packages.class_GRB import GRB
from packages.class_RSP import ResponseMatrix
from packages.package_simulations import simulate_observation

z = 1 # redshift 
imx, imy = 0, 0 # Position on the detector plane
ndets = 20000 # Number of enabled detectors

# Make a GRB object
template_grb = GRB()

# Make a Response Matrix object 
resp_mat = ResponseMatrix()

simulate_observation(template_grb, z, imx, imy, ndets, resp_mat)