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

		# Set light curve of GRB
		if light_curve_fn is not None:
			self.light_curve = fits.getdata(light_curve_fn,ext=1)

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