"""
Author: Mike Moss
Contact: mikejmoss3@gmail.com	

Test running sandbox and unit test runner

"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits 
import importlib
from util_packages.package_datatypes import dt_sim_res
from packages.package_bayesian_block import bayesian_t_blocks
from packages.class_SPECFUNC import PL, CPL
from packages.class_GRB import GRB
from packages.class_PLOTS import PLOTSIMRES, PLOTSAMPLE, PLOTGRB, PLOTS
from util_packages.package_det_ang_dependence import find_pcode, find_inc_ang, fraction_correction
from packages.package_many_simulations import many_simulations, make_param_list, make_ave_sim_res


grbs_names = np.array([
	"050416A",
	"050525A",
	"060614",
	"060912A",
	"061021",
	"080430",
	"080916A",
	"081007",
	"090424",
	"091018",
	# "091127",
	"100621A",
	"100625A",
	"100816A",
	"101219A",
	"110715A",
	"111228A",
	"120311A",
	"130427A",
	"130427A_cut",
	"130603B",
	"130925A",
	"140506A",
	"160425A",
	"160804A",
	"161001A",
	"161219B",
	], dtype="U11")




std1 = np.zeros(shape=0,dtype=[("name","U10"),("std",float)])
std64 = np.zeros(shape=0,dtype=[("name","U10"),("std",float)])

for i in range(len(grbs_names)):
	grbp = importlib.import_module("data_files.grb_{}.info".format(grbs_names[i]), package=None) # Load GRB parameters
	
	template_grb = GRB(grbname = grbp.name, z=grbp.z)
	template_grb.load_light_curve(grbp.fn, rm_trigtime=True, det_area=0.16)
	template_grb.cut_light_curve(tmin=-20, tmax=20)
	
	std = np.array( [(grbp.name, np.mean(template_grb.light_curve['UNC']) )] , dtype=[("name","U10"),("std",float)] )
	if template_grb.dt == 1:
		std1 = np.append(std1, std)
	if template_grb.dt < 1:
		std64 = np.append(std64, std)

print("For 1s ave = ", np.mean(std1['std']))
for i in range(len(std1)):
	print("\t{}, {:.2f}".format(std1['name'][i], std1['std'][i]))
print("For 64ms ave = ", np.mean(std64['std']))
for i in range(len(std64)):
	print("\t{}, {:.2f}".format(std64['name'][i], std64['std'][i]))


bins = np.linspace(start=np.min(std1['std']), stop=np.max(std1['std']), num=len(std1))
plt.hist(std1['std'],bins)
plt.show()