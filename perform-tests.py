"""
Author: Mike Moss
Contact: mikejmoss3@gmail.com	

Test running sandbox and unit test runner

"""

run_unit_tests = False 

if run_unit_tests is True:
	import unittest

	import unit_tests.test_class_SPECFUNC as test_sf
	import unit_tests.test_package_bayesian_block as test_bb

	runner = unittest.TextTestRunner()

	# Run unit tests:
	
	# Test classes
	runner.run(test_sf.suite())

	# Test packages
	runner.run(test_bb.suite())

	# Test utility packages

import numpy as np


a1 = np.ones(shape=(5,4))
a2 = (1,2,3,4)

print(a1*a2)