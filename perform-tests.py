"""
Author: Mike Moss
Contact: mikejmoss3@gmail.com	

Test running sandbox 

"""

# from unit_tests.test_bayesian_blocks import run_test
# run_test()

import numpy as np
import astropy
from astropy.modeling import Fittable1DModel, CompoundModel, Parameter

bins = 204 
emin = 15
emax = 350
spec_arr = np.zeros(shape=bins,dtype=[("ENERGY",float),("RATE",float),("UNC",float)])
spec_arr['ENERGY'] = np.logspace(np.log10(emin),np.log10(emax),num=bins)
model = astropy.modeling.powerlaws.PowerLaw1D(alpha=0.5,amplitude=1000)
spec_arr['RATE'] = model(spec_arr['ENERGY'])
spec_arr['UNC'] = np.sqrt(spec_arr['RATE'])
np.savetxt("data-files/test_spectrum.txt",spec_arr)