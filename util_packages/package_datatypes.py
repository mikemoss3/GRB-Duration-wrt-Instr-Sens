"""
Author: Mike Moss
Contact: mikejmoss3@gmail.com	

Defines data types commonly used in this code.
This makes it simple to change the shape of a structured array by only having it defined in one place.

"""

import numpy as np

# Data type for simulation results 
dt_sim_res = np.dtype([("DURATION",float),("TSTART",float),("FLUENCE",float),("z",float),("imx",float),("imy",float),("ndets",float)])