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
		self.cutoff_power_law = CPL()
		self.band = Band()

	def test_PL_eval(self):
		energy = 100 # energy to evaluate the power law at
		ans = 0.01 # expected flux photons / keV
		self.assertEqual(self.power_law(energy), ans, "Power law did not return correct flux.")

	def test_CPL_eval(self):
		energy = 100 # energy to evaluate the power law at
		ans = 0.01 # expected flux photons / keV
		self.assertEqual(self.cut_power_law(energy), ans, "Cut-off power law did not return correct flux.")

	def test_Band_eval(self):
		energy = 100 # energy to evaluate the power law at
		ans = 0.01 # expected flux photons / keV
		self.assertEqual(self.band(energy), ans, "Band function did not return correct flux.")

def suite():
	suite = unittest.TestSuite()
	suite.addTest(TestSpectralFunctions('test_PL_eval'))
	return suite