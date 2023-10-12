"""
Author: Mike Moss
Contact: mikejmoss3@gmail.com	

Test running sandbox and unit test runner

"""

run_unit_tests = True

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



# from packages.package_bayesian_block import bayesian_t_blocks
# from packages.class_GRB import GRB

# light_curve_fn = "./unit_tests/test_files/grb_130831A_1chan_64ms.txt"
# real_grb = GRB(light_curve_fn=light_curve_fn)

# duration, timestart, fluence = bayesian_t_blocks(real_grb,dur_per=90)
# print(duration,timestart,fluence)