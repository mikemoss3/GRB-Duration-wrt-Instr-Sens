# GRB Duration w.r.t. Instrument Sensitivity
Author: Mike Moss
Contact: mikejmoss3@gmail.com

## Purpose

This project allows a user to measure the duration of simulated gamma-ray burst (GRB) observations while taking into consideration observation conditions, such as the angle of the simulated GRB with respect to the detector bore-sight. This code is based on the work of [Moss et al. 2022](https://ui.adsabs.harvard.edu/abs/2022ApJ...927..157M/abstract).

## Procedure and How-to-Use

Here are short descriptions of each file and directory in the project:

* ** data-files/ **
	- Holds all input files to be used for and output files created from simulations
* ** packages/ **
	- Holds the core packages and classes needed to run the code
	__ class_GRB.py __
		- Defines a GRB object to store observed and simulated light curve and spectral information
	__ class_PLOTS.py __
		- Defines the class and methods used for plotting simulation results
	__ class_RSP.py __
		- Defines the main class this code uses to store response matrices and the associated methods
	__ class_SPECFUNC.py __
		- Defines the all classes related to spectral functions used in this library to represent GRB spectra
	__ package_analysis.py __
		- Defines functions to obtain the duration and fluence measurements for many synthetic GRBs
	__ package_bayesian_block.py __
		- Defines the Bayesian block method to calculate the duration of a GRB from a supplied light curve
	__ package_simulations.py __
		- Defines all the functions necessary to simulate an observation of a GRB using an input template, designated instrument response, and observing conditions
* ** unit_tests/ **
	- Holds all unit tests
* ** util_packages/ **
	- Holds the support packages and libraries for the main code 
* __ perform-analysis.py __
	- The intended main script of the code. All packages and classes can be called and used here.
* __ perform-tests.py __
	- A test-bed sandbox.


### Loading a Template GRB
First, create a GRB object that will act as a template for our simulations.
```
from packages.class_GRB import GRB

template_grb = GRB()
```

Next load a light curve for `template_grb`. Currently, the `GRB` class can load light curves from either .txt files or .fits files. In this example, the Swift/BAT light curve for GRB 081007 is used[^1]. The light curve is stored in the subdirectory `data-files/template-light-curves`
```
template_grb.load_light_curve("data-files/template-light-curves/grb_081007_1chan_1s.lc", rm_trigtime=True)
```

Now a spectral function must be defined for the GRB. In this example, a power law spectral function with a spectral index $\alpha = -1$ and $norm = 4$ (note: the normalization energy is set to $e_{norm} = 1$ keV by default[^2])
```
from packages.class_SPECFUNC import PL

spectral_function = PL(alpha=-1.,norm=4)
template_grb.load_specfunc( spectral_function )
```
This can be shortened to a single line like so 
```
template_grb.load_specfunc( PL(alpha=-1.,norm=4) )
```

Currently, the power law (PL), cut-off power law (CPL), and Band (Band) spectral functions are implemented.

[^1]: Light curves and spectral parameters for all Swift/BAT GRBs can be found on the online [Swift/BAT Catalog](https://swift.gsfc.nasa.gov/results/batgrbcat/)
[^2]: However, the spectral parameters found on the Swift/BAT catalog assume that the normalization energy is 50 keV (see the page 11 of the Third Swift/BAT GRB Catalog, [Lien et al. 2014](https://swift.gsfc.nasa.gov/results/batgrbcat/3rdBATcatalog.pdf))

### Loading an Instrument Response Matrix

### Simulating A GRB

### Applying Bayesian Blocks

### Simulating Many GRBs
Repeating the above steps for many observing condition combinations can be tedious, so the `package_analysis` package was developed to perform many simulations based on a given list of parameter combinations. Create a parameter list by defining the specific values of $z$, $imx$, $imy$, and $ndets$ desired,
```
from packages.package_analysis import make_param_list, many_simulations

z_arr = np.array([1])
imx_arr = np.array([0])
imy_arr = np.array([0])
ndets_arr = np.array([30000])
param_list = make_param_list(z_arr,imx_arr,imy_arr,ndets_arr)
```

Now, call the `many_simulations()` method. This requires specifying a template GRB that holds a user-defined light curve and spectral function, the parameter combination list that was just created, and a number of trials to simulate each parameter combination for. 
```
trials = 10
sim_results = many_simulations(template_grb, param_list, trials)
```
Three important keywords should be considered when running `many_simulations()` that are default to `False`. One, to run the code using multiple cores on your machine, set `multiproc=True`. This can greatly reduce the run time of the simulations. Two, if only the average duration of each parameter combination is wanted, then set `ret_ave=True`. Third, to keep a single light curve example for each parameter combination, set `keep_synth_grbs=True`. When setting `keep_synth_grbs` to `True`, a second output variable is required to hold the list of returned list of simulated GRB objects. All together this looks like
```
sim_results = many_simulations(template_grb, param_list, trials, multiproc=True, ret_ave=True, keep_synth_grbs=True)
```

### Plotting Simulation Results