"""
Author: Mike Moss
Contact: mikejmoss3@gmail.com	

Defines the unit tests for the Bayesian block package 

"""

from packages.package_bayesian_block import bayesian_t_blocks
from packages.class_GRB import GRB

def run_test():
	"""
	Run main test
	"""

	# Light curve file with known source
	light_curve_fn = "./unit_tests/test_files/grb_130831A_1chan_64ms.txt"
	# light_curve_fn = "./unit_tests/test_files/grb_160314A_1chan_64ms.txt"

	# Make GRB object to store light curve
	grb = GRB(light_curve_fn=light_curve_fn)

	# Run Bayesian block method
	duration, fluence = bayesian_t_blocks(grb)

	print(duration, fluence)

	# Run tests
	test_duration(duration)
	# test_fluence(fluence)

	return 0;

def test_duration(duration):
	"""
	This method test if the correct duration was recovered 
	"""
	known_duration = 8.732 # second

	if ((known_duration-0.5) < duration) and ( duration < (known_duration+0.5)):
		return 0
	else:
		print("Duration found by Bayesian blocks does not match known source duration.")
		return 1

def test_fluence(duration):
	"""
	This method test if the correct fluence was recovered 
	"""
	known_duration = 1 # second

	if (duration != known_duration):
		print("Fluence found by Bayesian blocks does not match known source duration.")
		return 1
	else:
		return 0
