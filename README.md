# GRB Duration w.r.t. Instrument Sensitivity
Author: Mike Moss
Contact: mikejmoss3@gmail.com

## Purpose

This project allows a user to measure the duration of simulated gamma-ray burst (GRB) observations while taking into consideration observation conditions, such as the angle of the simulated GRB with respect to the detector bore-sight. This code is based on the work of [Moss et al. 2022](https://ui.adsabs.harvard.edu/abs/2022ApJ...927..157M/abstract). This library currently focuses on analysis with the Swift/BAT instrument, but can easily be applied to any coded-aperture mask instrument. Support for non-coded-aperture mask instruments is in development! 

## Procedure and How-to-Use

Here are short descriptions of each file and directory in the project:
```
├── data-files/				# Holds all input files to be used for and output files created from simulations
├── packages/				# Holds the core packages and classes needed to run the code
│   ├── class_GRB.py				# Defines a GRB object to store observed and simulated light curve and spectral information
│   ├── class_PLOTS.py				# Defines the class and methods used for plotting simulation results
│   ├── class_RSP.py				# Defines the main class this code uses to store response matrices and the associated methods
│   ├── class_SPECFUNC.py			# Defines the all classes related to spectral functions used in this library to represent GRB spectra
│   ├── package_analysis.py			# Defines functions to obtain the duration and fluence measurements for many synthetic GRBs
│   ├── package_bayesian_block.py 		# Defines the Bayesian block method to calculate the duration of a GRB from a supplied light curve
│   ├── package_simulations.py			# Defines all the functions necessary to simulate an observation of a GRB using an input template, designated instrument response, and observing conditions
├── unit_tests/ 				# Holds all unit tests
├── util_packages/				# Holds the support packages and libraries for the main code 
├── perform-analysis.py				# The intended main script of the code. All packages and classes can be called and used here.
└── perform-tests.py				# A test-bed sandbox.
```

### Loading a Template GRB
First, create a GRB object that will act as a template for our simulations.
```
from packages.class_GRB import GRB

template_grb = GRB()
```

Next load a light curve for `template_grb`. Currently, the `GRB` class can load light curves from either .txt files or .fits files. In this example, the Swift/BAT mask-weighted light curve for GRB 081007 is used[^1]. The light curve is stored in the subdirectory `data-files/template-light-curves`
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

### Simulating A GRB
Before a synthetic GRB can be simulated, the desired observing condition parameter values should be assigned and a response matrix should be created. The following example will focus on the Swift/BAT instrument, but can rather easily be applied to any coded-aperture mask instrument.
```
from packages.class_RSP import ResponseMatrix

z_p = 1 # redshift to the synthetic GRB
imx = 0 # x-axis location of the source on the detector plane
imy = 0 # y-axis location of the source on the detector plane
ndets = 30000 # number of detectors enable on the detector plane at the time of observation

resp = ResponseMatrix() # Make response matrix object instance
```
There are many methods in place for loading response matrices. In the example below, a Swift/BAT response matrix will be loaded. The specific response matrix loaded depends on where GRB is located on the instrument detector plane (i.e., specified with the coordinates (imx, imy)). Instead of properly remaking the response matrix each time with the standard BAT tools, 30 response matrices have been stored in `/util_packages/files-swiftBAT-resp-mats/`. The response matrices roughly separate the BAT detector plane into 30 regions (see method from [Lien et al 2014](https://ui.adsabs.harvard.edu/abs/2014ApJ...783...24L/abstract)).
```
resp.load_SwiftBAT_resp(imx,imy)
```

To simulate a synthetic GRB observation and measure it's duration, import and run the command 
```
from packages.package_simulations import simulate_observation

synth_grb = simulate_observation( template_grb, z_p, imx, imy, ndets, resp)
```
And there you have it, the `synth_grb` object is a GRB object that has a mask-weighted light curve and spectrum generated from a known template input light curve and spectrum that have been folded through the Swift/BAT response. 

To obtain a duration measurement of the light curve, the `bayesian_t_blocks` method can be used. To find the T<sub>90</sub> duration, the duration percentage keyword should be set to 90, i.e., `dur_per=90` (note, this is the default value).
```
from packages.package_bayesian_block import bayesian_t_blocks

duration, t_start, fluence = bayesian_t_blocks(synth_grb, dur_per=90)
```

### Simulating Many GRBs
Repeating the above steps for many observing condition combinations can be tedious, so the `package_analysis` package was developed to perform many simulations based on a given list of parameter combinations. Create a parameter list by defining the specific values of $z$, $imx$, $imy$, and $ndets$ desired,
```
from packages.package_analysis import make_param_list, many_simulations

z_arr = np.array([1,2,3])
imx_arr = np.array([0,0.5])
imy_arr = np.array([0,0.5])
ndets_arr = np.array([30000,20000,10000])
param_list = make_param_list(z_arr,imx_arr,imy_arr,ndets_arr)
```
`param_list` now contains a list of every unique combination of parameters possible from the four parameter lists `z_arr`, `imx_arr`, `imy_arr`, and `ndets_arr`.

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