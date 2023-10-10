"""
Author: Mike Moss
Contact: mikejmoss3@gmail.com	

Defines the main class this code uses to store GRB observed and simulated data

"""

import numpy as np
from astropy.io import fits

from packages.class_SPECFUNC import SPECFUNC

class GRB(object):
	"""
	GRB class used to store observations information for observed and simulated GRBs

	Attributes:
	----------
	grbname : string
		Name of the GRB

	"""
	def __init__(self,grbname=None,z=0,imx=0,imy=0,
		T100_dur=None,T100_start=None,
		spectrum=None,light_curve_fn=None):

		# Assign this instance's parameters
		self.grbname = grbname
		self.z = z
		self.imx, self.imy = imx, imy
		self.T100_dur, self.T100_start = T100_dur, T100_start

		self.light_curve = None # Currently loaded light curve
		self.spectrum = SPECTRUM(model=None,params=None) # Currently loaded spectrum 
		self.spectra = np.zeros(shape=0,dtype=[("TSTART",float),("TEND",float),("SPECTRUM",SPECTRUM)]) # Time resolved spectrum array

		# Set light curve of GRB if a light curve file is supplied
		if light_curve_fn is not None:
			self.load_light_curve(light_curve_fn)
		# Set spectrum of GRB if a spectrum object is supplied
		if spectrum is not None:
			self.spectrum = spectrum

	def info_extractor(self,info_fn: str):
		"""
		Function to extract all information that is defined in a GRB information file.
		"""
		dtypes = np.dtype([('grb_name',"U100"),
			("lc","U100"),
			("z_0",float),
			("trials",int),
			("function","U100"),
			("alpha",float),
			("beta",float),
			("en_peak",float),
			("en_min",float),
			("en_max",float),
			("Quad","U10"),
			("pcode",float),
			("ndet",float),
			("rbsize",float),
			("spec_ev",int),
			("bckgrd","U100"),
			("det_trig_time",float),
			("T100_dur",float),
			("T100_start",float),
			("T100_end",float),
			("grb_flu",float),
			("theta",float),
			("imx",float),
			("imy",float),
			("GridID","U10")])
		info_arr = np.ndarray(shape=1,dtype=dtypes )

		with open(info_fn) as info:
			info_lines = info.readlines()
			for i in range(len(info_lines)):
				for j in range(len(dtypes)):
					if info_lines[i] ==  dtypes.names[j]+'=\n':
						# info_arr[0][j] = info_lines[i+1][0:-1]
						info_arr[dtypes.names[j]] = info_lines[i+1][0:-1]
		info.close()

		# We must ensure that alpha > -2, otherwise the Band function normalization will be equal to infinity.
		if float(info_arr['alpha']) < -2.0:
			info_arr['alpha'] = -1.99

	def load_spectrum(self,model,params,tstart=None,tend=None):
		"""
		Method to load a spectrum

		Attributes:
		----------
		model : string
			Name of the spectral model to load
		params : list
			List of the spectral parameters for this model
		tstart, tend : float, float
			Used to indicate the start and stop time of a time-resolved spectrum. If None are given, a time-average spectrum is assumed.
		"""

		# Time resolved spectrum
		if tstart is not None:
			# Check that both a start and stop time were given 
			if tend is None:
				print("Please provide both a start and end time.")
				return 0;

			# Check if this is the first loaded spectrum 
			if len(self.spectra) == 0:
				self.spectra = np.insert(self.spectra,0,(tstart,tend,SPECTRUM(model,params)))
				return 0;
			else:
				# If not, find the index where to insert this spectrum (according to the time)
				for i in range(len(self.spectra)):
					if self.spectra[i]['TSTART'] > tstart:
						# Insert the new spectrum 
						self.spectra = np.insert(self.spectra,i,(tstart,tend,SPECTRUM(model,params)))
						return 0;
					# If the new spectrum is the last to start, append it to the end
					self.spectra = np.insert(self.spectra,len(self.spectra),(tstart,tend,SPECTRUM(model,params)))
					return 0;
		# Time averaged spectrum
		else:
			self.spectrum = SPECTRUM(model,params)
			
			return 0;

	def make_spectrum(self,emin,emax,num_bins = None,spec_num=None):
		"""
		Method to evaluate the spectrum over the defined energy interval using the GRB object's spectral model and parameters

		Attributes:
		----------
		emin, emax : float, float
			Defines the lower and upper bounds of the energy interval over which to evaluate the spectrum
		num_bins : int
			Number of energy bins to use, default is 10*log(emax/emin)
		"""

		if num_bins is None:
			num_bins = np.log10(emax/emin)*10

		if spec_num is None:
			model = self.spectrum.model
			params = self.spectrum.params
		else:
			model = self.spectra[spec_num].model
			params = self.spectra[spec_num].params

		# Initialize array
		spectrum = np.zeros(shape=num_bins,dtype=[("ENERGY",float),("RATE",float)])
		# Evaluate energy array
		spectrum['ENERGY'] = np.logspace(np.log10(emin),np.logspace(emax),num=num_bins)

		# Evaluate spectrum rate
		spectrum['RATE'] = model(spectrum['ENERGY'],*params)

		return spectrum

	def load_light_curve(self,file_name,inc_unc = True,t_offset=0,rm_trigtime=False):
		
		# Check if this is a fits file or a text file 

		if file_name.endswith(".lc") or file_name.endswith(".fits"):
			if inc_unc is False:
				tmp_light_curve = fits.getdata(file_name,ext=1)
				self.light_curve = np.zeros(shape=len(tmp_light_curve),dtype=[('TIME',float),('RATE',float)])
				self.light_curve['TIME'] = tmp_light_curve['TIME']
				if rm_trigtime is True:
						self.light_curve['TIME']-=fits.getheader(file_name,ext=0)['TRIGTIME']
				self.light_curve['RATE'] = tmp_light_curve['RATE']
			else:
				tmp_light_curve = fits.getdata(file_name,ext=1)
				self.light_curve = np.zeros(shape=len(tmp_light_curve),dtype=[('TIME',float),('RATE',float),('UNC',float)])
				self.light_curve['TIME'] = tmp_light_curve['TIME']
				if rm_trigtime is True:
						self.light_curve['TIME']-=fits.getheader(file_name,ext=0)['TRIGTIME']
				self.light_curve['RATE'] = tmp_light_curve['RATE']
				self.light_curve['UNC'] = tmp_light_curve['ERROR']
		elif file_name.endswith(".txt"):
			if inc_unc is False:
				self.light_curve = np.genfromtxt(file_name,dtype=[('TIME',float),('RATE',float)])
			else:
				self.light_curve = np.genfromtxt(file_name,dtype=[('TIME',float),('RATE',float),('UNC',float)])

		if t_offset != 0:
			self.light_curve -= t_offset

	def move_to_source_frame(self,z, rm_bgd_sig=True):
		"""
		Method to place the GRB light curve and spectra into the source rest frame using the given redshift

		Attributes:
		----------
		z : float
			Redshift to correct for 
		"""

		# Check if a spectrum is loaded
		if (self.spectrum._check_model()==1) or (self.spectrum._check_params()==1):
			return;

		# Remove background signal outside of T100 (requires that the T100 start time and duration were defined)
		if rm_bgd_sig is True:
			inds = np.where( (self.light_curve['TIME'] < self.T100_start) & (self.light_curve['TIME'] > (self.T100_start+self.T100_dur)) )

		## 
		# Shift Spectrum
		## 
		


		##
		# Shift Light Curve
		## 

		# Apply distance corrections to flux values
		
		# Apply time-dilation to light curve (i.e., correct time-binning)


