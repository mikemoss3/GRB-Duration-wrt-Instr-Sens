"""
Author: Mike Moss
Contact: mikejmoss3@gmail.com	

Test running sandbox and unit test runner

"""

from astropy.io import fits 
from packages.package_bayesian_block import custom_bblocks, bayesian_t_blocks
from packages.class_GRB import GRB

grb= GRB()
# grb.load_light_curve("data_files/grb_060614/grb_060614_1chan_1s.lc", rm_trigtime=True, det_area=0.16)
grb.load_light_curve("data_files/grb_211211A/grb_211211A_1chan_1s.lc", rm_trigtime=True, det_area=0.16)
grb.cut_light_curve(tmin=-100, tmax=250)

# print(bayesian_t_blocks(grb))

duration, fluence = custom_bblocks(grb.light_curve)
print(duration)


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
