import numpy as np
import subprocess
from subprocess import STDOUT

obs_low_z_grbs = np.array([
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
	"091127",
	"100621A",
	"100625A",
	"100816A",
	"101219A",
	"110715A",
	"111228A",
	"120311A",
	"130427A",
	"130603B",
	"130925A",
	"140506A",
	"160425A",
	"160804A",
	"161001A",
	"161219B",
	], dtype="U10")

for i in range(len(obs_low_z_grbs)):
	subprocess.run(["cp ../grb_{}/grb_{}_redshift_sim-results.tmp.txt.npy ./grb_{}_redshift_sim-results.txt.npy".format(obs_low_z_grbs[i], obs_low_z_grbs[i], obs_low_z_grbs[i])], shell=True, stderr=STDOUT)
	subprocess.run(["cp ../grb_{}/grb_{}_redshift_ave-sim-results.tmp.txt.npy ./grb_{}_redshift_ave-sim-results.txt.npy".format(obs_low_z_grbs[i], obs_low_z_grbs[i], obs_low_z_grbs[i])], shell=True, stderr=STDOUT)