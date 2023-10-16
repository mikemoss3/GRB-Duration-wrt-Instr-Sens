# GRB Duration w.r.t. Instrument Sensitivity
Author: Mike Moss
Contact: mikejmoss3@gmail.com

## Purpose

This project allows a user to quickly simulate the observation of gamma-ray bursts (GRBs) while taking into consideration observation conditions, such as the angle of the simulated GRB with respect to the detector bore-sight.

## Procedure and How-to-Use

First, create a GRB object that will act as a template for our simulations.
```
from packages.class_GRB import GRB

template_grb = GRB()
```
Then load a light curve for `template_grb`. Currently, the `GRB` class can load light curves from either .txt files or .fits files. In this example, the Swift/BAT light curve for GRB 081007 is used[^1].
```
template_grb.load_light_curve("data-files/template-light-curves/grb_081007_1chan_1s.lc", rm_trigtime=True)
```
Now a spectral function must be defined for the GRB. In this example, a power law spectral function with a spectral index $\alpha = -1$ and $norm = 4$ (note: the normalization is set to $e_{norm} = 1$ keV by default[^2])
```
from packages.class_SPECFUNC import PL

spectral_function = PL(alpha=-1.,norm=4)
template_grb.load_specfunc( spectral_function )
```
or 
```
template_grb.load_specfunc( PL(alpha=-1.,norm=4) )
```
Currently, the power law (PL), cut-off power law (CPL), and Band (Band) spectral functions are implemented.

[^1]: Light curves and spectral parameters for all Swift/BAT GRBs can be found on the online [Swift/BAT Catalog](https://swift.gsfc.nasa.gov/results/batgrbcat/)
[^2]: However, the spectral parameters found on the Swift/BAT catalog assume that the normalization energy is 50 keV (see the page 11 of the Third Swift/BAT GRB Catalog, [Lien et al. 2014](https://swift.gsfc.nasa.gov/results/batgrbcat/3rdBATcatalog.pdf))