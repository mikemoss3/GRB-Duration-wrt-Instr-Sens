import numpy as np
from scipy.integrate import romberg

def make_spectrum(model,model_params,emin,emax,en_bins):
	"""
	Make a spectrum from E_min to E_max using the specified model
	"""


	spectrum = np.zeros(shape=en_bins,dtype=[("ENERGY",float),("RATE",float)])
	spectrum['ENERGY'] = np.logspace(np.log(emin),np.log(emax),en_bins)
	spectrum['RATE'] = model(spectrum['ENERGY'],*model_params)

	return spectrum

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

def _energy_flux(Emin, Emax,specmodel,params):
	"""
	Method to find the total energy flux of the spectral model within the given energy range

	Input Emin and Emax in units of keV
	"""
	energy_flux_kev = romberg(function=specmodel,args=params,a=Emin,b=Emax)  ## [keV/s/cm2]

	kev2erg = 1000*1.60217657e-12

	energy_flux = energy_flux_kev*kev2erg  ## [erg/s/cm2]

	return energy_flux

def _find_norm(flux,Emin,Emax,specmodel,params):
	"""
	Method to find the spectrum normalization based on observed flux
	"""
	return flux/_energy_flux(Emin,Emax,specmodel,params)