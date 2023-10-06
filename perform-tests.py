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

alpha = 1.0
norm = 2.06201**(-2)

spec_arr = np.zeros(shape=bins,dtype=[("ENERGY",float),("RATE",float),("UNC",float)])
spec_arr['ENERGY'] = np.logspace(np.log10(emin),np.log10(emax),num=bins)
model = astropy.modeling.powerlaws.PowerLaw1D(alpha=alpha,amplitude=2000)
spec_arr['RATE'] = model(spec_arr['ENERGY']) * norm
spec_arr['UNC'] = np.sqrt(spec_arr['RATE']) * norm
np.savetxt("data-files/test_spectrum.txt",spec_arr)

# from packages_util.package_det_ang_dependence import find_pcode
# print(find_pcode(-0.24,0.))