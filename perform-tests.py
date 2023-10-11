"""
Author: Mike Moss
Contact: mikejmoss3@gmail.com	

Test running sandbox 

"""

# from unit_tests.test_bayesian_blocks import run_test
# run_test()

import numpy as np
import matplotlib.pyplot as plt

from packages.class_GRB import GRB
from packages.class_SPECFUNC import PL

bins = 204 
emin = 15
emax = 350

alpha = -1.0
norm = 2.06201**(-2)

# Make a GRB object
template_grb = GRB()
# Make light curve 
template_grb.load_light_curve("data-files/sw00330856000b_1chan_64ms.lc", rm_trigtime=True)
template_grb.light_curve = template_grb.light_curve[np.argmax(-100 <= template_grb.light_curve['TIME']):np.argmax(template_grb.light_curve['TIME'] >= 100)]

template_grb.load_spectrum(PL(alpha=alpha,norm=norm))


plt.plot(template_grb.light_curve['TIME'],template_grb.light_curve['RATE'])
template_grb.move_to_new_frame(0.5, 0.2,rm_bgd_sig=False)
plt.plot(template_grb.light_curve['TIME'],template_grb.light_curve['RATE'])


plt.show()
