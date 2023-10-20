"""
Author: Mike Moss
Contact: mikejmoss3@gmail.com	

Test running sandbox and unit test runner

"""

import matplotlib.pyplot as plt 
import numpy as np

from packages.class_GRB	import GRB
from packages.class_SPECFUNC import PL, CPL
from packages.class_RSP import ResponseMatrix
from packages.class_PLOTS import PLOTS

from util_packages.package_det_ang_dependence import fraction_correction

plot = PLOTS()

template_grb = GRB(z=0.1)
template_grb.load_light_curve("data-files/template-light-curves/grb_160121A_1chan_1s.lc", rm_trigtime=True) # counts / sec / det
template_grb.light_curve = template_grb.light_curve[np.argmax(-100 <= template_grb.light_curve['TIME']):np.argmax(template_grb.light_curve['TIME'] >= 100)]
template_grb.light_curve['RATE'] /= 0.16  # counts / sec / cm^2
template_grb.light_curve['UNC'] /= 0.16

alpha = -1.77
# tmin = -11.156
# tmax = 36.76
# norm = template_grb.get_ave_photon_flux(tmin=tmin,tmax=tmax) / PL(alpha= alpha,norm=1,enorm=1)._calc_phot_flux(15, 150) # for GRB GRB 081007
# norm = 0.6 / PL(alpha= alpha,norm=1,enorm=1)._calc_phot_flux(15, 150) # for GRB GRB 081007
# print(template_grb.get_ave_photon_flux(tmin=tmin,tmax=tmax))
# print(norm * (50)**alpha)
# template_grb.load_specfunc(PL(alpha= alpha,norm=norm,enorm=1))
norm = 4.432*10**(-3)
frac = fraction_correction(0.417,-0.056)
template_grb.load_specfunc(PL(alpha= alpha,norm=norm/frac/11,enorm=50))

# ax = plt.figure().gca()

resp = ResponseMatrix()
resp.load_rsp_from_file('util_packages/files-swiftBAT-resp-mats/grb_160121A_preslew.rsp')
# resp.plot_effarea()
# resp.load_SwiftBAT_resp(0.,0.)
# resp.plot_effarea(ax=ax)

ax_spec = plt.figure().gca()
plot.plot_spectra(template_grb,emin=10,emax=150,en_window=[10,150],ax=ax_spec)
plot.plot_spectra(template_grb,resp,en_window=[10,150],ax=ax_spec)

# plt.savefig("data-files/figs/2023-10-20/grb-160121A-spectra.png")

# template_grb.light_curve['RATE']*=0.16
# template_grb.light_curve['UNC']*=0.16
# ax_l_c = plt.figure().gca()
# plot.plot_light_curves(template_grb,ax=ax_l_c)

# plt.savefig("data-files/figs/2023-10-20/grb-160121A-light-curve.png")


plt.show()



run_unit_tests = False
if run_unit_tests is True:
	import unittest

	import unit_tests.test_class_SPECFUNC as test_sf

	import unit_tests.test_package_bayesian_block as test_bb
	import unit_tests.test_package_cosmology as test_cos
	import unit_tests.test_package_det_ang_dependence as test_dad

	runner = unittest.TextTestRunner()

	# Run unit tests:
	
	# Test classes
	runner.run(test_sf.suite())

	# Test packages
	runner.run(test_bb.suite())
	runner.run(test_cos.suite())
	runner.run(test_dad.suite())


