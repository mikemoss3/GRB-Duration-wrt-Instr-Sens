"""
Author: Mike Moss
Contact: mikejmoss3@gmail.com	

Defines the main class this code uses to store response matrices and the associated methods.
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.stats import norm
import copy 

from util_packages.package_det_ang_dependence import find_grid_id


class ResponseMatrix(object):
	"""
	Response matrix class
	"""
	def __init__(self,E_phot_min=1,E_phot_max=500,num_phot_bins=200,E_chan_min=1,E_chan_max=200,num_chans=80):
		"""
		RSP class initialization 
		Attributes:
		E_phot_min = 0, minimum photon energy, keV 
		E_phot_max = 500, maximum photon energy, keV
		num_phot_bins = 200, number of bins along the photon energy axis
		E_chan_min = 0, minimum channel energy, keV
		E_chan_max = 200, maximum channel energy, keV
		num_chans = 80, number of instrument channels 
		"""

		self.num_phot_bins = num_phot_bins
		self.num_chans = num_chans
		
		self.set_E_phot(E_phot_min=E_phot_min, E_phot_max=E_phot_max, num_phot_bins=num_phot_bins,verbose=False)
		self.set_E_chans(E_chan_min=E_chan_min,E_chan_max=E_chan_max,num_chans=num_chans,verbose=False)
		
		self.N_GRP = np.ones(shape=num_phot_bins) # The number of 'channel subsets' for for the energy bin
		self.F_CHAN = np.zeros(shape=num_phot_bins) # The channel number of the start of each "channel subset" for the energy bin
		self.N_CHAN = np.ones(shape=num_phot_bins)*num_chans # The number of channels within each "channel subset" for the energy bin

		# Initialize self.MATRIX as empty array
		self.make_empty_resp() # Contains all the response probability values for each
		# 						'channel subset' corresponding to the energy bin for a given row


	def __copy__(self):
		cls = self.__class__
		result = cls.__new__(cls)
		result.__dict__.update(self.__dict__)
		return result

	def __deepcopy__(self, memo):
		cls = self.__class__
		result = cls.__new__(cls)
		memo[id(self)] = result
		for k, v in self.__dict__.items():
			setattr(result, k, copy.deepcopy(v, memo))
		return result

	def copy(self):
		return copy.deepcopy(self)

	def deepcopy(self):
		return copy.deepcopy(self)

	def set_E_phot(self,E_phot_min=0,E_phot_max=500,num_phot_bins=200,verbose=True):
		""" Class method to set the photon energy axis""" 
		ENERG_AX = make_en_axis(E_phot_min, E_phot_max, num_phot_bins)
		self.num_phot_bins = num_phot_bins
		self.ENERG_LO = ENERG_AX['Elo'] # Incoming photon energy, lower bound
		self.ENERG_MID =  ENERG_AX['Emid'] # Incoming photon energy, bin center
		self.ENERG_HI =  ENERG_AX['Ehi'] # Incoming photon energy, upper bound
		if verbose is True:
			print("Response matrix has been reset to zeros.")
		self.make_empty_resp()

	def set_E_chans(self,E_chan_min=0,E_chan_max=200,num_chans=80,verbose=True):
		""" Class method to set the photon energy axis""" 
		ECHAN_AX = make_en_axis(E_chan_min, E_chan_max, num_chans)
		self.num_chans = num_chans
		self.ECHAN_LO = ECHAN_AX['Elo'] # Instrument energy channel lower bound
		self.ECHAN_MID = ECHAN_AX['Emid'] # Instrument energy channel center
		self.ECHAN_HI = ECHAN_AX['Ehi'] # Instrument energy channel upper bound
		if verbose is True:
			print("Response matrix has been reset to zeros.")
		self.make_empty_resp()

	def make_empty_resp(self):
		""" If the shape of the response matrix is changed, the response matrix is reset to zeros."""
		self.MATRIX = np.zeros(shape=(self.num_phot_bins,self.num_chans))

	def identity(self):
		""" Make identity matrix """
		for i in range(self.num_phot_bins):
			for j in range(self.num_chans):
				if i==j:
					self.MATRIX[i,j] = 1

	def overDeltaE(self,alpha=2):
		""" Decrease as 1/DeltaE^alpha from E_true """
		for i in range(self.num_phot_bins):
			for j in range(self.num_chans):
				self.MATRIX[i,j] = 1/(1+np.abs(self.ECHAN_MID[j] - self.ENERG_MID[i])**alpha)
			
			# Normalize this column
			self.MATRIX[i]/=np.sum(self.MATRIX[i])

	def gauss(self):
		""" Decrease as probability in a Gaussian behavior as channel energy moves away from photon energy """
		for i in range(self.num_phot_bins):
			dist = norm(self.ENERG_MID[i],np.sqrt(self.ENERG_MID[i]))
			for j in range(self.num_chans):
				self.MATRIX[i,j] = dist.pdf(self.ECHAN_MID[j])

	def load_rsp_from_file(self,file_name):
		""" Load response matrix from file""" 
		resp_data = fits.getdata(filename=file_name,extname="SPECRESP MATRIX")
		ebounds_data = fits.getdata(file_name,extname="EBOUNDS")
		
		self.num_phot_bins = len(resp_data)
		self.num_chans = len(ebounds_data)

		self.ENERG_LO = np.zeros(shape=self.num_phot_bins)
		self.ENERG_HI = np.zeros(shape=self.num_phot_bins)
		self.ENERG_MID = np.zeros(shape=self.num_phot_bins)
		self.N_GRP = np.zeros(shape=self.num_phot_bins)
		self.F_CHAN = np.zeros(shape=self.num_phot_bins)
		self.N_CHAN = np.zeros(shape=self.num_phot_bins)
		self.MATRIX = np.zeros(shape=(self.num_phot_bins,self.num_chans) )

		self.ECHAN_LO = np.zeros(shape=self.num_chans)
		self.ECHAN_HI = np.zeros(shape=self.num_chans)
		self.ECHAN_MID = np.zeros(shape=self.num_chans)

		for i in range(self.num_phot_bins):
			self.ENERG_LO[i] = resp_data[i][0] # Incoming photon energy, lower bound
			self.ENERG_HI[i] =  resp_data[i][1] # Incoming photon energy, upper bound
			self.ENERG_MID[i] =  (self.ENERG_LO[i]+self.ENERG_HI[i])/2 # Incoming photon energy, bin center
			self.N_GRP[i] = resp_data[i][2] # The number of 'channel subsets' for for the energy bin
			self.F_CHAN[i] = resp_data[i][3] # The channel number of the start of each "channel subset" for the energy bin
			self.N_CHAN[i] = resp_data[i][4] # The number of channels within each "channel subset" for the energy bin
			self.MATRIX[i] = resp_data[i][5] # Contains all the response probability values for each
			# 										'channel subset' corresponding to the energy bin for a given row
		
		for i in range(self.num_chans):
			self.ECHAN_LO[i] = ebounds_data[i][1] # Instrument energy channel lower bound
			self.ECHAN_HI[i] = ebounds_data[i][2] # Instrument energy channel upper bound
			self.ECHAN_MID[i] = (self.ECHAN_LO[i]+self.ECHAN_HI[i])/2 # Instrument energy channel center

	def load_SwiftBAT_resp(self,imx,imy):
		"""
		Method to load an (averaged) Swift/BAT response matrix given the position of the source on the detector plane.
		"""

		# Obtain GridID
		gridid = find_grid_id(imx,imy)

		# Load corresponding response matrix
		self.load_rsp_from_file(file_name = "./util_packages/files-swiftBAT-resp-mats/BAT_alldet_grid_{}.rsp".format(gridid))


	def plot_heatmap(self,ax=None,E_phot_bounds=None,E_chan_bounds=None):
		""" Plot heat map of the response matrix """

		if ax is None:
			ax = plt.figure().gca()
		fig = plt.gcf()

		im = ax.pcolormesh(self.ECHAN_MID,self.ENERG_MID,self.MATRIX,shading='auto')

		if E_chan_bounds is None:
			ax.set_xlim(self.ECHAN_HI[0],self.ECHAN_LO[-1])
		else:
			ax.set_xlim(E_chan_bounds[0],E_chan_bounds[1])
		if E_phot_bounds is None:
			ax.set_ylim(self.ENERG_HI[0],self.ENERG_LO[-5])
		else:
			ax.set_xlim(E_phot_bounds[0],E_phot_bounds[1])

		ax.set_xlabel('Instrument Channel Energy (keV)')
		ax.set_ylabel('Photon Energy (keV)')

		cbar = fig.colorbar(im)
		cbar.ax.set_ylabel('Probability', rotation=270,labelpad=15)

	def plot_effarea(self,ax=None,det_area=1,E_phot_bounds=None,norm=1):
		""" Plot heat map of the response matrix """
		
		if ax is None:
			ax = plt.figure().gca()

		# eff_area = np.sum(self.MATRIX,axis=1)/(self.ENERG_HI-self.ENERG_LO)
		eff_area = np.zeros(shape=len(self.MATRIX))
		for i in range(len(self.MATRIX)):
			for j in range(len(self.MATRIX[0])):
				eff_area[i] += self.MATRIX[i][j]
		
		eff_area*=det_area

		ax.step(self.ENERG_MID,eff_area*norm)

		if E_phot_bounds is None:
			ax.set_xlim(self.ENERG_MID[0],self.ENERG_MID[-1])
		else:
			ax.set_xlim(E_phot_bounds[0],E_phot_bounds[1])

		ax.set_xscale('log')
		# ax.set_yscale('log')

		ax.set_xlabel('Incident Photon Energy (keV)')
		ax.set_ylabel(r'Effective Area (cm$^2$)')

	def fold_spec(self,specfunc):
		"""
		Method to fold a spectrum through this response matrix
		"""

		folded_spec = make_folded_spec(specfunc,self)

		return folded_spec


# class ResponseMatrixArray(ResponseMatrix):
# 	"""
# 	Class to store response matrices for more than one time interval
# 	"""
# 	def __init__(self):

def make_en_axis(Emin,Emax,num_en_bins):
	""" Make energy axis """

	en_axis = np.zeros(shape=num_en_bins,dtype=[("Elo",float),("Emid",float),("Ehi",float)])
	en_axis['Elo'] = np.logspace(np.log10(Emin),np.log10(Emax),num_en_bins,endpoint=False)
	en_axis['Ehi'] = np.logspace(np.log10(en_axis["Elo"][1]),np.log10(Emax),num_en_bins,endpoint=True)
	en_axis['Emid'] = (en_axis['Ehi'] + en_axis['Elo'])/2
	return en_axis

def make_folded_spec(source_spec_func,rsp):
	""" 
	Convolve spectral function with instrument response to obtain observed spectrum

	Attributes:
	----------
	source_spec_func : SPECFUNC
		Unfolded source spectral function
	rsp : RSP
		Response matrix 
	"""

	# Initialize folded spectrum 
	folded_spec = np.zeros(shape=rsp.num_chans,dtype=[("ENERGY",float),("RATE",float),("UNC",float)])
	# The folded spectrum will have the same energy bins as the response matrix
	folded_spec['ENERGY'] = rsp.ECHAN_MID

	# Initialize the binned source spectrum
	# If the source spectrum covers a smaller interval than the response matrix, any energy bin outside of the source spectrum energy range will have a rate equal to zero.
	binned_source_spec = np.zeros(shape=rsp.num_phot_bins,dtype=[("ENERGY",float),("RATE",float)])

	binned_source_spec['ENERGY'] = rsp.ENERG_MID
	binned_source_spec['RATE'] = source_spec_func(binned_source_spec['ENERGY'])

	# Fold the correctly binned source spectrum with the response matrix
	folded_spec['RATE'] = np.matmul(binned_source_spec['RATE'],rsp.MATRIX)/(rsp.ECHAN_HI - rsp.ECHAN_LO)

	# What should the uncertainty be?
	# folded_spec['UNC'] = np.sqrt(folded_spec['RATE'])
	folded_spec['UNC'] = 0.05*folded_spec['RATE']

	return folded_spec
