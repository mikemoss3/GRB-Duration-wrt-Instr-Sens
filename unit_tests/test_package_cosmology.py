import unittest

from util_packages.package_cosmology import lum_dis, k_corr
from packages.clss_SPECFUNC import PL, CPL

class TestCosmology(unittest.TestCase):
	def test_lum_dis(self):
		z = 1  # input redshift
		known_dis = 1  # cm, known luminosity distance
		
		self.assertEqual(lum_dis(z), known_dis, "Luminosity distance is incorrect.")

	def test_k_corr_PL(self):
		"""
		Test the k-correction for a power law spectral function
		"""
		power_law = PL(alpha=-1,norm=1)  # Instance of a power law spectral function 
		z = 1  # Redsihft
		emin = 15  # Energy band minimum
		emax = 350  # Energy band maximum 

		known_k_corr = 1.  # Known k-correction value

		self.assertEqual(k_corr(power_law,z,emin,emax), known_k_corr,"k-correction of a power law is incorrect.")

	def test_k_corr_CPL(self):
		"""
		Test the k-correction for a cut-off power law spectral function
		"""
		cutoff_power_law = CPL(alpha=-1,ep=300,norm=1)  # Instance of a cut-off power law spectral function 
		z = 1  # Redsihft
		emin = 15  # Energy band minimum
		emax = 350  # Energy band maximum 

		known_k_corr = 1.  # Known k-correction value

		self.assertEqual(k_corr(cutoff_power_law,z,emin,emax), known_k_corr,"k-correction of a cut-off power law is incorrect.")


def suite():
	suite = unittest.TestSuite()
	suite.addTest(TestCosmology('test_lum_dis'))
	return suite