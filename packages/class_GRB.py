"""
Author: Mike Moss
Contact: mikejmoss3@gmail.com	

Defines the main class this code uses to store GRB observed and simulated data

"""

import numpy as np
from astropy.io import fits
import copy 

from packages.class_SPECFUNC import SPECFUNC
from util_packages.package_cosmology import lum_dis, k_corr
import util_packages.globalconstants as gc


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
		self.specfunc = None # Currently loaded spectrum 
		self.spectrafuncs = np.zeros(shape=0,dtype=[("TSTART",float),("TEND",float),("SPECFUNC",SPECFUNC)]) # Time resolved spectrum array

		# Set light curve of GRB if a light curve file is supplied
		if light_curve_fn is not None:
			self.load_light_curve(light_curve_fn)
		# Set spectrum of GRB if a spectrum object is supplied
		if spectrum is not None:
			self.spectrum = spectrum

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

	def load_specfunc(self,specfunc,tstart=None,tend=None):
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
			if len(self.spectrafuncs) == 0:
				self.spectrafuncs = np.insert(self.spectrafuncs,0,(tstart,tend,specfunc))
				return 0;
			else:
				# If not, find the index where to insert this spectrum (according to the time)
				for i in range(len(self.spectrafuncs)):
					if self.spectrafuncs[i]['TSTART'] > tstart:
						# Insert the new spectrum 
						self.spectrafuncs = np.insert(self.spectra,i,(tstart,tend,specfunc))
						return 0;
					# If the new spectrum is the last to start, append it to the end
					self.spectrafuncs = np.insert(self.spectra,len(self.spectra),(tstart,tend,specfunc))
					return 0;
		# Time averaged spectrum
		else:
			self.specfunc = specfunc
			
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
			specfunc = self.specfunc
		else:
			specfunc = self.specfunc[spec_num]

		spectrum = specfunc.make_spectrum(emin,emax,num_bins)

		return spectrum

	def load_light_curve(self,file_name,inc_unc = True,t_offset=0,rm_trigtime=False,T100_dur=None,T100_start=None):
		"""
		Method to load a light curve from either a .fits or .txt file
		"""

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
			self.light_curve['TIME'] -= t_offset
		if T100_dur is not None:
			self.T100_dur = T100_dur
		if T100_start is not None:
			self.T100_start = T100_start

	def move_to_new_frame(self, z_o, z_p, emin=gc.bol_lum[0],emax=gc.bol_lum[1],rm_bgd_sig=False):
		"""
		Method to shift the GRB light curve and spectra from a frame at z_o to a frame at z_p

		if z_p = 0, this is the same as shifting the GRB to the source frame and the light curve returned will be the bolometric one.
		If z_o = 0, it is assumed the GRB is in the rest frame

		Attributes:
		----------
		z_o : float
			Current redshift of the GRB
		z_p : float
			Redshift to shift the GRB to
		emin, max : float, float
			Spectrum energy band minimum and maximum
		rm_bgd_sig : bool
			Indicates whether or not to remove the background signal outside the T100 range should be removed. 
		"""
		if z_o == z_p:
			# No chage in the light curve or spectrum.
			return;

		# Update redshift class attribute
		self.z = z_p

		# Remove background signal outside of T100 (requires that the T100 start time and duration were defined)
		if rm_bgd_sig is True:
			inds = np.where( (self.light_curve['TIME'] > self.T100_start) & (self.light_curve['TIME'] < (self.T100_start+self.T100_dur)) )
			new_light_curve = np.zeros(shape=len(self.light_curve))
			new_light_curve[inds] = self.light_curve[inds]

		# Calculate k_correction factor
		kcorr = k_corr(self.specfunc, z_o, emin, emax)/k_corr(self.specfunc, z_p, emin, emax)
		# Move spectral function to z_p frame by correcting E_peak or temperature by the redshift (if spectral function has a peak energy or temperature)
		for i, (key, val) in enumerate(self.specfunc.params.items()):
			if key == "ep":
				self.specfunc.params[key] *= (1+z_o)/(1+z_p)
			if key == "temp":
				self.specfunc.params[key] *= (1+z_o)/(1+z_p)

		##
		# Shift Light Curve
		## 

		# Apply distance corrections to flux values (See Bloom, Frail, and Sari 2001 Equation 4)
		dis_corr_to_z_o, dis_corr_to_z_p = 1., 1.
		if z_o != 0:
			dis_corr_to_z_o = 4 * np.pi * np.power(lum_dis(z_o), 2.) / (1+z_o)
		if z_p != 0:
			dis_corr_to_z_p = 4 * np.pi * np.power(lum_dis(z_p), 2.) / (1+z_p)
		self.light_curve['RATE'] = self.light_curve['RATE'] * kcorr * dis_corr_to_z_o / dis_corr_to_z_p
		self.light_curve['UNC'] = self.light_curve['UNC'] * kcorr * dis_corr_to_z_o / dis_corr_to_z_p
		
		# Apply time-dilation to light curve (i.e., correct the time binning)
		# Calculate the start and stop times of the flux light curve in the z_p frame.
		tpstart = self.light_curve['TIME'][0]*(1+z_p)/(1+z_o)
		tpend = self.light_curve['TIME'][-1]*(1+z_p)/(1+z_o)

		# Bin size of the light curve curve
		bin_size = (self.light_curve['TIME'][1] - self.light_curve['TIME'][0])
		# Create a time axis from tpstart to tpend with proper bin size
		tmp_time_arr = np.arange(tpstart,tpend+bin_size,bin_size)

		# Create an array to store the flux light curve in the z_p frame
		flux_lc_at_z_p = np.zeros(shape=len(tmp_time_arr),dtype=([("TIME",float),("RATE",float)]))
		flux_lc_at_z_p['TIME'] = tmp_time_arr

		# Temporary light curve to store z_p frame light curve
		tmp_light_curve = np.zeros(shape=len(tmp_time_arr),dtype = [("TIME",float),("RATE",float),("UNC",float)])
		tmp_light_curve['TIME'] = tmp_time_arr

		# We must correct for time dilation by binning the flux into z_p frame time bins
		if z_o > z_p:
			# The light curve must be squeezed

			# For each time bin of the z_p light curve curve:
			for i in range(len(tmp_light_curve)-1):
				# What time bin is this?
				curr_time_bin = tmp_light_curve['TIME'][i]

				# Find the indices of the z_o frame light curve where the time axis is encompassed by curr_time_bin*(1+z_o)/(1+z_p) and (1+curr_time_bin)*(1+z_o)/(1+z_p)
				argstart = np.argmax(self.light_curve['TIME']>=curr_time_bin*(1+z_o)/(1+z_p))
				argend = np.argmax(self.light_curve['TIME']>=(curr_time_bin+bin_size)*(1+z_o)/(1+z_p))

				# Grab the total sum of the distance corrected flux within these time bins.
				tmp_light_curve['RATE'][i] = np.sum(self.light_curve['RATE'][argstart:argend])
				tmp_light_curve['UNC'][i] = np.sum(np.power(self.light_curve['UNC'][argstart:argend],2)) # Add uncertainties in quadrature 
		else:
			# z_p > z_o, the light curve must be stretched

			# For each time bin of the z_o light curve:
			for i in range(len(self.light_curve)-1):
				# What time bin is this?
				curr_time_bin = self.light_curve['TIME'][i]

				# Grab the distance corrected flux from the z_o. This flux must be distributed across multiple time bins.
				curr_flux_to_distribute = self.light_curve['RATE'][i]
				curr_flux_unc_to_distribute = self.light_curve['UNC'][i]

				# Find the indices of z_p light curve where the time axis encompasses curr_time_bin*(1+z_p)/(1+z_o) and (1+curr_time_bin)*(1+z_p)/(1+z_o)
				argstart = np.argmax(tmp_light_curve['TIME']>=curr_time_bin*(1+z_p)/(1+z_o))
				argend = np.argmax(tmp_light_curve['TIME']>=(curr_time_bin+bin_size)*(1+z_p)/(1+z_o))

				# How many time bins are within the time interval of curr_time_bin*(1+z_p)/(1+z_o) and (1+curr_time_bin)*(1+z_p)/(1+z_o) ?
				# We need to know the number of time bins because we must divide the curr_flux_to_distribute across each time bin in the z_p frame
				# that came from the (curr_time_bin : 1+curr_time_bin) time interval in the z_o frame.
				num_new_time_bins = len(tmp_light_curve['TIME'][argstart:argend])
				if num_new_time_bins == 0:
					argend = argstart+1
					num_new_time_bins = 1

				# For the bins between (argstart:argend) assign the flux value of curr_flux_to_distribute/num_new_time_bins
				tmp_light_curve['RATE'][argstart:argend] = tmp_light_curve['RATE'][argstart:argend] + np.ones(shape=num_new_time_bins)*(curr_flux_to_distribute/num_new_time_bins)
				tmp_light_curve['UNC'][argstart:argend] = tmp_light_curve['UNC'][argstart:argend] + np.ones(shape=num_new_time_bins)*(curr_flux_unc_to_distribute/num_new_time_bins)

		# Set the light curve to the distance corrected light curve
		self.light_curve = tmp_light_curve

		return;
