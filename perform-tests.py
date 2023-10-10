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

from packages.class_SPECTRUM import SPECTRUM
from packages.class_SPECTRUM import PL, CPL, Band

test_spec = SPECTRUM(model=PL,params=[-alpha,2000*norm])

test_spec_arr = test_spec.make_spectrum(emin=emin,emax=emax,num_bins=bins)

import matplotlib.pyplot as plt
plt.plot(test_spec_arr['ENERGY'],test_spec_arr['RATE'])

plt.xscale('log')
plt.yscale('log')
plt.show()
