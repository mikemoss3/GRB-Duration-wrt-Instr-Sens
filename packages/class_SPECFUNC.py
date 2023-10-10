import numpy as np
from scipy.integrate import romberg

class SPECFUNC():
	"""
	Base class to hold variables associated with the loaded spectral model and its parameters

	Attributes
	----------
	model : string
		Name of the spectral model to use. Current available models are Band (Band), Power Law (PL), or Cut-off Power Law (CPL)
	params : list
		List of model parameters
	"""
	def __init__(self,params=None,**kwargs):
		
		self.color = "k"  # Define a color, used for plotting
		self.param_names = list(self.param_dict.keys())

		# If no parameter value list is given, use defaults.
		if params is None:
			self.params = list(self.param_dict.values())
		# Else, update the parameter values
		else: 
			self.params = params
			for i, key in self.param_names.keys():
				self.param_dict[key] = params[i]

		# If specific keyword parameters are defined
		for i, (key, val) in enumerate(kwargs.items()):
			self.param_dict[key] = arg

	def __call__(self,energy):
		"""
		Method to evaluate the spectral function at a given energy.
		"""
		return self.evaluate(energy,*self.params)

	def print_params(self):
		"""
		Print current parameters
		"""

		for param in self.param_names:
			print("Name: {}\n\tdescription: {}\n\tvalue = {}\n".format(param,getattr(type(self), param)._description,getattr(type(self), param).value))

		return 0;

	def model_name(self):
		"""
		Name of the compound model, automatically set as the combination of submodel names
		"""
		return self.name;

class PL(SPECFUNC):
	"""
	Power Law

	Parameters
	----------
	alpha : float
		Power law index
	norm : float
		Model normalization
	enorm : float
		Normalization energy
	"""
	name = "Power Law"

	alpha = -1.
	norm = 1.
	enorm = 1.

	param_dict = {"alpha" : alpha ,"norm" : norm ,"enorm" : enorm}

	@staticmethod
	def evaluate(energy, alpha, norm, enorm=1):
		"""
		Compute the power law spectrum at a particular energy given the current spectral parameters
		"""

		flux_value = norm * np.power(energy/enorm, alpha)
		
		return flux_value

	def testmeth(self):
		print(self.alpha)

def CPL(SPECFUNC):
	"""
	Cut-off Power Law 

	Parameters
	----------
	ep : float
		Peak energy
	alpha : float
		Power law index
	norm : float
		Model normalization
	enorm : float
		Normalization energy
	"""
	name = "Cut-off Power Law"
	param_dict = {"ep","alpha","norm","enorm = 1 keV"}

	ep = 100.
	alpha = -1.
	norm = 1.
	enorm = 1.

	@staticmethod
	def evaluate(energy, ep, alpha, norm, enorm):
		"""
		Compute the cut-off power law spectrum at a particular energy given the current spectral parameters
		"""
		flux_value = norm * np.power(energy/enorm, alpha) * np.exp(- energy / ep)

		return flux_value

class Blackbody(SPECFUNC):
	"""
	Blackbody function.

	Parameters
	----------
	temp : float
		Blackbody temperature (in units of energy, i.e., k_B*T where k_B is the Boltzmann constant)
	alpha : float
		Index of the power law below temperature
	norm : float
		Model normalization
	"""
	name = "Blackbody"

	temp = 20.
	alpha = 0.4
	norm = 1

	@staticmethod
	def evaluate(energy, temp, alpha, norm):
		"""
		Compute the blackbody spectrum at a particular energy
		"""

		# Initialize the return value
		flux_value = np.zeros_like(energy,subok=False)

		i = energy < 2e3
		if i.max():
			# If the energy is less than 2 MeV
			flux_value[i] = norm * np.power(energy[i]/temp,1.+alpha)/(np.exp(energy[i]/temp) - 1.)
		i = energy >= 2e3
		if i.max():
			flux_value[i] = 0

		return flux_value

def Band(SPECFUNC):
	"""
	Band function (see Band et al. 1993)

	Parameters
	----------
	ep : float
		Peak energy
	alpha : float
		Low energy power law index
	beta : float
		High energy power law index
	norm : float
		Model normalization
	"""
	name = "Band"
	param_dict = {"ep","alpha","beta","norm","enorm = 100 keV"}

	ep = 100.
	alpha = -1.
	beta = -2.
	norm = 1.
	enorm = 1.

	@staticmethod
	def evaluate(energy, ep, alpha, beta, norm,enorm=100):
		"""
		Compute the Band spectrum at a particular energy given the current spectral parameters
		"""
		# Initialize the return value
		flux_value = np.zeros_like(energy,subok=False)

		# Calculate break energy
		e0 = ep / (alpha - beta)

		test

		i = energy <= ep
		if i.max():
			flux_value[i] = norm * np.power(energy[i]/enorm, alpha) * np.exp(- energy[i] / e0)
		
		i = energy > ep
		if i.max():
				flux_value[i] = norm * np.power((alpha - beta) * e0/enorm, alpha - beta) * np.exp(beta - alpha) * np.power(energy[i]/enorm,beta)

		return flux_value




# """
# def make_spectrum(self,emin,emax,num_bins = None):
# 		"""
# 		Method to evaluate the spectrum over the defined energy interval using the GRB object's spectral model and parameters

# 		Attributes:
# 		----------
# 		emin, emax : float, float
# 			Defines the lower and upper bounds of the energy interval over which to evaluate the spectrum
# 		num_bins : int
# 			Number of energy bins to use, default is 10*log(emax/emin)
# 		"""
# 		# Check if model and parameters are loaded
# 		if (self._check_model()==1) or (self._check_params()==1):
# 			return;

# 		if num_bins is None:
# 			num_bins = np.log10(emax/emin)*10

# 		# Initialize array
# 		spectrum = np.zeros(shape=num_bins,dtype=[("ENERGY",float),("RATE",float)])
# 		# Evaluate energy array
# 		spectrum['ENERGY'] = np.logspace(np.log10(emin),np.log10(emax),num=num_bins)

# 		# Evaluate spectrum rate
# 		spectrum['RATE'] = self.model(spectrum['ENERGY'],*self.params)

# 		return spectrum

# 	def _energy_flux(self,emin,emax):
# 		"""
# 		Method to find the total energy flux of the spectral model within the given energy range

# 		Input Emin and Emax in units of keV
# 		"""
# 		energy_flux_kev = romberg(function=self.model,args=self.params,a=emin,b=emax)  ## [keV/s/cm2]

# 		kev2erg = 1000*1.60217657e-12

# 		energy_flux = energy_flux_kev*kev2erg  ## [erg/s/cm2]

# 		return energy_flux

# 	def _find_norm(self,flux,emin,emax):
# 		"""
# 		Method to find the spectrum normalization based on observed flux
# 		"""
# 		return flux/self._energy_flux(emin,emax)

# """