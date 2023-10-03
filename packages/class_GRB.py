"""
Author: Mike Moss
Contact: mikejmoss3@gmail.com	

Defines the main class this code uses to store GRB observed and simulated data

"""

import numpy as np
from astropy.io import fits

class GRB(object):
	"""
	GRB class used to store observations information for observed and simulated GRBs

	Attributes:

	"""
	def __init__(self,grbname=None,z=0,imx=0,imy=0,
		spec_func=None,alpha=None,beta=None,en_peak=None,
		T100_dur=None,T100_start=None,light_curve_fn=None):

		# Assign this instance's parameters
		self.grbname = grbname
		self.z = z
		self.imx, self.imy = imx, imy
		self.spec_func = spec_func
		self.alpha, self.beta, self.en_peak = alpha, beta, en_peak
		self.T100_dur, self.T100_start = T100_dur,T100_start

		self.light_curve = None # Currently loaded light curve
		self.spectrum = None # Currently loaded spectrum 
		self.spectra = np.zeros(shape=0,dtype=[("TSTART",float),("TEND",float),("SPECTRUM",tuple)]) # Time resolved spectrum array

		# Set light curve of GRB
		if light_curve_fn is not None:
			self.load_light_curve(light_curve_fn)

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

	def load_spectrum(self,file_name,tstart=None,tend=None):

		if tstart is not None:
			# Check that both a start and stop time were given 
			if tend is None:
				print("Please provide both a start and end time (or neither).")
				return 0;

			# Load the spectrum designated by the file name
			
			loaded_spectrum = np.genfromtxt(file_name,dtype=[("ENERGY",float),("RATE",float),("UNC",float)])

			# Check if this is the first loaded spectrum 
			if len(self.spectra) == 0:
				self.spectra = np.insert(self.spectra,0,(tstart,tend,loaded_spectrum))
				return 0;
			else:
				# If not, find the index where to insert this spectrum (according to the time)
				for i in range(len(self.spectra)):
					if self.spectra[i]['TSTART'] > tstart:
						# Insert the new spectrum 
						self.spectra = np.insert(self.spectra,i,(tstart,tend,loaded_spectrum))
						return 0;
					# If the new spectrum is the last to start, append it to the end
					self.spectra = np.insert(self.spectra,len(self.spectra),(tstart,tend,loaded_spectrum))
					return 0;
		else:
			self.spectrum = np.genfromtxt(file_name,dtype=[("ENERGY",float),("RATE",float),("UNC",float)])
			return 0;

	def set_spectrum(self,spec_array,tstart=None,tend=None):

		# Check that the array is in the correct format
		if spec_array.dtype.names != ('ENERGY','RATE','UNC'):
			print("Please provide an array with three columns: ENERGY, RATE, and UNC (all floats)")

		if tstart is not None:
			# Check that both a start and stop time were given 
			if tend is None:
				print("Please provide both a start and end time.")
				return 0;

			# Check if this is the first loaded spectrum 
			if len(self.spectra) == 0:
				self.spectra = np.insert(self.spectra,0,(tstart,tend,spec_array))
				return 0;
			else:
				# If not, find the index where to insert this spectrum (according to the time)
				for i in range(len(self.spectra)):
					if self.spectra[i]['TSTART'] > tstart:
						# Insert the new spectrum 
						self.spectra = np.insert(self.spectra,i,(tstart,tend,spec_array))
						return 0;
					# If the new spectrum is the last to start, append it to the end
					self.spectra = np.insert(self.spectra,len(self.spectra),(tstart,tend,spec_array))
					return 0;
		else:
			self.spectrum = spec_array
			return 0;

	def make_spectrum(self,emin,emax,bins,model,params,spec_num=None):
		"""
		Method to make the GRB spectrum using a supplied spectral model and its parameters
		"""

		# Initialize array
		if spec_num == None:
			self.spectrum = np.zeros(shape=bins,dtype=[("ENERGY",float),("RATE",float)])
		# If spec_num is not None:

		self.spectrum['ENERGY'] = np.logspace(np.log10(emin),np.logspace(emax),num=bins)
		self.spectrum['RATE'] = model(*params)


	def load_light_curve(self,file_name,inc_unc = True,t_offset=0):
		
		# Check if this is a fits file or a text file 

		if file_name.endswith(".lc") or file_name.endswith(".fits"):
			if inc_unc is False:
				tmp_light_curve = fits.getdata(file_name,ext=1)
				self.light_curve = np.zeros(shape=len(tmp_light_curve),dtype=[('TIME',float),('RATE',float)])
				self.light_curve['TIME'] = tmp_light_curve['TIME']
				self.light_curve['RATE'] = tmp_light_curve['RATE']
			else:
				tmp_light_curve = fits.getdata(file_name,ext=1)
				self.light_curve = np.zeros(shape=len(tmp_light_curve),dtype=[('TIME',float),('RATE',float),('UNC',float)])
				self.light_curve['TIME'] = tmp_light_curve['TIME']
				self.light_curve['RATE'] = tmp_light_curve['RATE']
				self.light_curve['UNC'] = tmp_light_curve['ERROR']
		elif file_name.endswith(".txt"):
			if inc_unc is False:
				self.light_curve = np.genfromtxt(file_name,dtype=[('TIME',float),('RATE',float)])
			else:
				self.light_curve = np.genfromtxt(file_name,dtype=[('TIME',float),('RATE',float),('UNC',float)])

		if t_offset != 0:
			self.light_curve -= t_offset
