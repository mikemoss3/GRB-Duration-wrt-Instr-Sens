"""
Author: Mike Moss
Contact: mikejmoss3@gmail.com	

Test running sandbox and unit test runner

"""

import matplotlib.pyplot as plt 
import numpy as np

from packages.class_GRB	import GRB
from packages.class_SPECFUNC import PL, CPL
from packages.class_RSP import ResponseMatrix


grb = GRB(z=0.1)
grb.load_light_curve("data-files/template-light-curves/grb_160121A_1chan_1s.lc", rm_trigtime=True) # counts / sec / det
grb.light_curve = grb.light_curve[np.argmax(-100 <= grb.light_curve['TIME']):np.argmax(grb.light_curve['TIME'] >= 100)]
grb.light_curve['RATE'] /= 0.16  # counts / sec / cm^2
grb.light_curve['UNC'] /= 0.16

grb.load_specfunc(PL(alpha=-1.77,norm=1,enorm=50)) # Photons / sec / keV / cm^2


ax = plt.figure().gca()

resp = ResponseMatrix()
resp.load_rsp_from_file('sw00671231000b_preslew.rsp')
resp.plot_effarea(ax=ax,norm = 15580)
resp.load_SwiftBAT_resp(0.,0.)
resp.plot_effarea(ax=ax)
plt.show()

# folded_spec = resp.fold_spec(spec)



run_unit_tests = False
if run_unit_tests is True:
	import unittest

	import unit_tests.test_class_SPECFUNC as test_sf

	import unit_tests.test_package_bayesian_block as test_bb
	import unit_tests.test_package_cosmology as test_cos
	import unit_tests.test_package_det_ang_dependence as test_dad

	runner = unittest.TextTestRunner()

	# Run unit tests:
	
	# Test classes
	runner.run(test_sf.suite())

	# Test packages
	runner.run(test_bb.suite())
	runner.run(test_cos.suite())
	runner.run(test_dad.suite())


