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
from packages.class_PLOTS import PLOTS

from util_packages.package_det_ang_dependence import fraction_correction, find_grid_id, find_pcode
from packages.package_many_simulations import many_simulations, make_param_list, make_ave_sim_res


rsp = ResponseMatrix()
rsp.load_SwiftBAT_resp(0, 0)
# rsp.load_SwiftBAT_resp(1.2, .75)

rsp.plot_heatmap()
# rsp.plot_effarea()

plt.show()

"""
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
"""
