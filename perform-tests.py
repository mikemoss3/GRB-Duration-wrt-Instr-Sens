"""
Author: Mike Moss
Contact: mikejmoss3@gmail.com	

Test running sandbox and unit test runner

"""


# import numpy as np
# from packages.class_SPECFUNC import PL, CPL, Band

# spec = Band()
# print(spec(5))

# energy = 10 
# energy = np.array([10,20,30,40,50,60])
# e_peak = 50
# i = energy < e_peak

# print(i)
# print(energy[i])



run_unit_tests = True
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


