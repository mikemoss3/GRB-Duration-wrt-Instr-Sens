import numpy as np

from packages.package_bayesian_block import bayesian_blocks
from packages.class_GRB import GRB

def run_test():
	"""
	Run main test
	"""

	# Light curve file with known source
	light_curve_fn = "light_curve.txt"

	# Make GRB object to store light curve
	grb = GRB(light_curve_fn=light_curve_fn)

	# Run Bayesian block method
	duration, fluence = bayesian_blocks(grb)

	# Run tests
	test_duration(duration)
	test_fluence(fluence)

	return 0;

def test_duration(duration):
	"""
	This method test if the correct duration was recovered 
	"""
	known_duration = 1 # second

	if (duration != known_duration):
		print("Duration found by Bayesian blocks does not match known source duration.")
		return 1
	else:
		return 0

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


if __name__ == "__main__":
	run_test()