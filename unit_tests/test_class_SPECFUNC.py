"""
Author: Mike Moss
Contact: mikejmoss3@gmail.com	

Defines the unit tests for the SPECFUNC class

"""

import unittest

from packages.class_SPECFUNC import PL, CPL, Band


class TestSpectralFunctions(unittest.TestCase):
	def setUp(self):
		self.power_law = PL()

	def test_PL_eval(self):
		energy = 100 # energy to evaluate the power law at
		ans = 0.01 # expected flux photons / keV
		self.assertEqual(self.power_law(energy), ans, "Power law did not return correct energy.")

	# test powerlaw
	# test cutoffpowerlaw
	# test band

def suite():
	suite = unittest.TestSuite()
	suite.addTest(TestSpectralFunctions('test_PL_eval'))
	return suite