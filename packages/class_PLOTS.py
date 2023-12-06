"""
Author: Mike Moss
Contact: mikejmoss3@gmail.com	

Defines the class and methods used for plotting simulation results.

"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

class PLOTS(object):
	"""
	Class that defines all methods used for plotting simulation results

	Attributes:
	----------
	sim_results : np.ndarray
		Array of simulation results
	"""

	def __init__(self,sim_results=None,
		fontsize = 13, fontweight = "bold"):
		self.sim_results = sim_results

		self.fontsize = fontsize
		self.fontweight = fontweight

	def plot_aesthetics(self,ax,xax=True,yax=True):
		"""
		This function is used to make bold and increase the font size of all plot tick markers

		Attributes:
		----------
		ax : matplotlib.axes
			Axis object to apply aesthetic changes to
		xax : boolean
			Indicates whether to apply to changes to x-axis
		yax : boolean
			Indicates whether to apply to changes to y-axis
		"""

		fontsize = self.fontsize
		fontweight = self.fontweight

		if xax is True:
			for tick in ax.xaxis.get_major_ticks():
				tick.label1.set_fontsize(fontsize=fontsize)
				tick.label1.set_fontweight(fontweight)

				tick.label2.set_fontsize(fontsize=fontsize)
				tick.label2.set_fontweight(fontweight)

			for tick in ax.xaxis.get_minor_ticks():
				tick.label1.set_fontweight(fontweight)

				tick.label2.set_fontweight(fontweight)

		if yax is True:
			for tick in ax.yaxis.get_major_ticks():
				tick.label1.set_fontsize(fontsize=fontsize)
				tick.label1.set_fontweight(fontweight)

				tick.label2.set_fontsize(fontsize=fontsize)
				tick.label2.set_fontweight(fontweight)

			for tick in ax.yaxis.get_minor_ticks():
				tick.label1.set_fontweight(fontweight)

				tick.label2.set_fontweight(fontweight)
			
		ax.tick_params(direction="in",which="both")
		ax.margins(x=0,y=0)

	def duration_overlay(self,light_curve,order_type=2,ax=None,**kwargs):
		"""
		Method to plot simulated duration measures overlayed on template light curve

		Attributes:
		----------
		light_curve : np.ndarray
			Light curve to plot under the simulated duration measurements
		ax : matplotlib.axes
			Axis on which to create the figure
		"""

		if ax is None:
			ax = plt.figure().gca()

		# Plot light curve
		ax.step(light_curve['TIME'],light_curve['RATE'],color="k",alpha=0.5,**kwargs)

		# Order sim_results
		# Take lines and order them
		# 0 == No order
		# 1 == Time Start
		# 2 == Time Duration
		if order_type == 0:
			sorted_sim_results = self.sim_results
		elif order_type == 1:
			sorted_sim_results = np.sort(self.sim_results,order='TSTART')
		elif order_type == 2:
			sorted_sim_results = np.flip(np.sort(self.sim_results,order="DURATION"))

		# Plot simulated duration measurements
		y_pos = np.linspace(np.max(light_curve['RATE'])*0.05,np.max(light_curve['RATE'])*0.95,len(sorted_sim_results))
		for i in range(len(sorted_sim_results)):
			ax.hlines(y=y_pos[i],xmin=sorted_sim_results[i]['TSTART'],xmax=(sorted_sim_results[i]['TSTART']+sorted_sim_results[i]['DURATION']),color="C1",alpha=0.7)

		ax.set_xlabel("Time (sec)",fontsize=self.fontsize,fontweight=self.fontweight)
		ax.set_ylabel("Rate (counts/sec)",fontsize=self.fontsize,fontweight=self.fontweight)

		self.plot_aesthetics(ax)

	def dur_vs_param(self,obs_param,t_true=None,dur_frac=False,ax=None,marker=".",**kwargs):
		"""
		Method to plot duration vs observing parameter (e.g., redshift, pcode, ndets)

		Attributes:
		----------
		obs_param : str
			The sim_result column field name of the parameter to be plotted against the duration (e.g., "z", "imx", or "ndets")
		t_true : float
			Value of the true duration of the burst. If given, a horizontal line will be marked at t_true.
		dur_frac : boolean
			Indicates whether the y-axis will be simply the duration measure (dur_frac = False) or the fraction of the true duration (dur_frac = True). If True, t_true must be supplied
		ax : matplotlib.axes
			Axis on which to create the figure
		"""

		if dur_frac is True:
			if t_true is None:
				print("A true duration must be given to create duration fraction axis.")
				return;
			sim_results = self.sim_results
			sim_results['DURATION'] /= t_true
		else:
			sim_results = self.sim_results

		if ax is None:
			ax = plt.figure().gca()

		ax.scatter(sim_results[obs_param],sim_results['DURATION'],marker=marker,**kwargs)

		if (t_true is not None):
			if (dur_frac is False):
				ax.axhline(y=t_true,color="C2",linestyle="dashed",alpha=0.5,label="True Duration")
			else:
				ax.axhline(y=1,color="C2",linestyle="dashed",alpha=0.5,label="True Duration")

		ax.set_xlabel("{}".format(obs_param),fontsize=self.fontsize,fontweight=self.fontweight)
		ax.set_ylabel("Duration (sec)",fontsize=self.fontsize,fontweight=self.fontweight)

		self.plot_aesthetics(ax)
		ax.margins(x=0.1,y=0.05)

	def det_plane_map(self,ax=None,imx_max=1.5,imy_max=1,dimx=0.1,dimy=0.1,**kwargs):
		"""
		Method to plot the average duration percentage as a function of the position on the detector plane

		Attributes:
		----------
		ax : matplotlib.axes
			Axis on which to create the figure
		imx_max, imy_max : float, float
			Defines the maximum (and minimum) values of the x and y plane on the detector
		"""

		if ax is None:
			ax = plt.figure().gca()
		fig = plt.gcf()

		im = ax.scatter(self.sim_results['imx'],self.sim_results['imy'],c=self.sim_results['DURATION'],cmap='viridis',**kwargs)
		cbar = fig.colorbar(im)

		ax.axhline(y=0,color="k",alpha=0.2)
		ax.axvline(x=0,color="k",alpha=0.2)

		ax.set_xlim(-imx_max,imx_max)
		ax.set_ylim(-imy_max,imy_max)

		ax.set_xlabel("IMX",fontsize=self.fontsize,fontweight=self.fontweight)
		ax.set_ylabel("IMY",fontsize=self.fontsize,fontweight=self.fontweight)

		cbar.set_label("Duration (sec)",fontsize=self.fontsize,fontweight=self.fontweight)

		fig.tight_layout()
		self.plot_aesthetics(ax)

	def plot_light_curves(self,grbs,t_window=None,labels=None,ax=None,alpha=0.7,**kwargs):
		"""
		Method to plot the average duration percentage as a function of the position on the detector plane

		Attributes:
		----------
		grbs : GRB, array of GRB
			Either a single instance of a GRB or an array of GRBs for which the light curves will be plotted on the same axis
		ax : matplotlib.axes
			Axis on which to create the figure
		"""

		if ax is None:
			ax = plt.figure().gca()

		# For an array of GRBs
		if hasattr(grbs,'__len__'):
			for i in range(len(grbs)):
				if labels is None:
					ax.errorbar(x=grbs[i].light_curve['TIME'],y=grbs[i].light_curve['RATE']*grbs[i].dt,yerr=grbs[i].light_curve['UNC']*grbs[i].dt,fmt="",drawstyle="steps-mid",alpha=alpha,**kwargs)
				else:
					ax.errorbar(x=grbs[i].light_curve['TIME'],y=grbs[i].light_curve['RATE']*grbs[i].dt,yerr=grbs[i].light_curve['UNC']*grbs[i].dt,fmt="",drawstyle="steps-mid",alpha=alpha,label="{}".format(labels[i]),**kwargs)
		# For a single GRB
		else:
			ax.errorbar(x=grbs.light_curve['TIME'],y=grbs.light_curve['RATE']*grbs.dt,yerr=grbs.light_curve['UNC']*grbs.dt,fmt="",drawstyle="steps-mid",alpha=alpha,label=labels,**kwargs)

		ax.set_xlabel("Time (sec)",fontsize=self.fontsize,fontweight=self.fontweight)
		ax.set_ylabel("Rate (counts/sec)",fontsize=self.fontsize,fontweight=self.fontweight)

		if t_window is not None:
			ax.set_xlim(t_window)

		if labels is not None:
			ax.legend(fontsize=self.fontsize-2)

		self.plot_aesthetics(ax)


	def plot_spectra(self,grbs,resp=None,emin=None,emax=None,en_window=None,labels=None,ax=None,alpha=0.7,norm=1,**kwargs):
		"""
		Method to plot the average duration percentage as a function of the position on the detector plane

		Attributes:
		----------
		grbs : GRB, array of GRB
			Either a single instance of a GRB or an array of GRBs for which the spectra will be plotted on the same axis
		resp : RSP
			Response matrix object to fold the source spectra with in order to produce a folded spectrum.
		emin, emax : float, float
			Minimum and maximum energy to plot the source spectrum over.
		ax : matplotlib.axes
			Axis on which to create the figure
		"""


		if ax is None:
			ax = plt.figure().gca()

		# For an array of GRBs
		if hasattr(grbs,'__len__'):
			for i in range(len(grbs)):
				if resp is None:
					spectrum = grbs[i].make_spectrum(emin,emax)
					ax.step(x=spectrum['ENERGY'],y=spectrum['RATE']*norm,alpha=alpha,label="{}".format(labels[i]),**kwargs)
				else: 
					folded_spec = resp.fold_spec(grbs[i].specfunc) 
					ax.step(x=folded_spec['ENERGY'],y=folded_spec['RATE']*norm,alpha=alpha,label="{}".format(labels[i]),**kwargs)
		# For a single GRB
		else:
			if resp is None:
				spectrum = grbs.make_spectrum(emin,emax)
				ax.step(x=spectrum['ENERGY'],y=spectrum['RATE']*norm,alpha=alpha,label=labels,**kwargs)
			else: 
				folded_spec = resp.fold_spec(grbs.specfunc) 
				ax.step(x=folded_spec['ENERGY'],y=folded_spec['RATE']*norm,alpha=alpha,label=labels,**kwargs)

		ax.set_xscale('log')
		ax.set_yscale('log')

		if en_window is not None:
			ax.set_xlim(en_window)
	
		ax.set_xlabel("Time (sec)",fontsize=self.fontsize,fontweight=self.fontweight)
		ax.set_ylabel("Rate (counts/sec)",fontsize=self.fontsize,fontweight=self.fontweight)

		if labels is not None:
			ax.legend(fontsize=self.fontsize-2)

		self.plot_aesthetics(ax)