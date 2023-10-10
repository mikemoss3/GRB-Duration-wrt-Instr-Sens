import numpy as np
from scipy.integrate import romberg

class SPECTRUM(object):
	"""
	Class to hold variables associated with the loaded spectral model and its parameters

	Attributes
	----------
	model : string
		Name of the spectral model to use. Current available models are Band (Band), Power Law (PL), or Cut-off Power Law (CPL)
	params : list
		List of model parameters
	"""
	def __init__(self,model,params):
		self.model = model
		self.params = params

	def make_spectrum(self,emin,emax,num_bins = None):
		"""
		Method to evaluate the spectrum over the defined energy interval using the GRB object's spectral model and parameters

		Attributes:
		----------
		emin, emax : float, float
			Defines the lower and upper bounds of the energy interval over which to evaluate the spectrum
		num_bins : int
			Number of energy bins to use, default is 10*log(emax/emin)
		"""
		# Check if model and parameters are loaded
		if (self._check_model()==1) or (self._check_params()==1):
			return;

		if num_bins is None:
			num_bins = np.log10(emax/emin)*10

		# Initialize array
		spectrum = np.zeros(shape=num_bins,dtype=[("ENERGY",float),("RATE",float)])
		# Evaluate energy array
		spectrum['ENERGY'] = np.logspace(np.log10(emin),np.log10(emax),num=num_bins)

		# Evaluate spectrum rate
		spectrum['RATE'] = self.model(spectrum['ENERGY'],*self.params)

		return spectrum

	def _energy_flux(self,emin,emax):
		"""
		Method to find the total energy flux of the spectral model within the given energy range

		Input Emin and Emax in units of keV
		"""
		energy_flux_kev = romberg(function=self.model,args=self.params,a=emin,b=emax)  ## [keV/s/cm2]

		kev2erg = 1000*1.60217657e-12

		energy_flux = energy_flux_kev*kev2erg  ## [erg/s/cm2]

		return energy_flux

	def _find_norm(self,flux,emin,emax):
		"""
		Method to find the spectrum normalization based on observed flux
		"""
		return flux/self._energy_flux(emin,emax)

	def _check_model(self):
		"""
		Test if a model has been loaded
		"""
		if self.model is None:
			print("Error: A model has not been loaded.")
			return 1;
		return 0;

	def _check_params(self):
		"""
		Test if parameters have been loaded
		"""
		if self.params is None:
			print("Error: Model parameters has not been loaded.")
			return 1;
		return 0;

def Band(energy, e0, alpha, beta, norm,enorm=100):
	"""
	Compute the Band function spectrum at a particular energy or array of energies

	Parameters
	----------
	e0 : float
		Break energy
	alpha : float
		Low energy power law index
	beta : float
		High energy power law index
	norm : float
		Model normalization
	"""

	# Initialize the return value
	flux_value = np.zeros_like(energy,subok=False)

	# Calculate peak energy
	e_lim = (alpha - beta)*e0

	i = energy <= e_lim
	if i.max():
		flux_value[i] = norm * np.power(energy[i]/enorm, alpha) * np.exp(- energy[i] / e0)
	
	i = energy > e_lim
	if i.max():
			flux_value[i] = norm * np.power((alpha - beta) * e0/enorm, alpha - beta) * np.exp(beta - alpha) * np.power(energy[i]/enorm,beta)

	return flux_value

def PL(energy, alpha, norm, enorm=1):
	"""
	Compute the power law function spectrum at a particular energy or array of energies

	Parameters
	----------
	alpha : float
		Power law index
	norm : float
		Model normalization
	enorm : float
		Normalization energy
	"""

	flux_value = norm * np.power(energy/enorm, alpha)
	
	return flux_value

def CPL(energy, eb, alpha, norm, enorm=1):
	"""
	Compute the power law function spectrum at a particular energy or array of energies

	Parameters
	----------
	eb : float
		Break energy
	alpha : float
		Power law index
	norm : float
		Model normalization
	enorm : float
		Normalization energy
	"""

	flux_value = norm * np.power(energy/enorm, alpha) * np.exp(- energy / enorm)

	return flux_value