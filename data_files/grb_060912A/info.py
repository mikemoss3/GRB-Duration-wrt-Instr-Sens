name = "060912A"
z = 0.94 # Measured redshift 
zmax = 5.5 # Highest redshift to simulate at
fn = "data_files/grb_060912A/grb_060912A_1chan_64ms.lc" # file path to light curve
t_true = 11.5 # true T90
t_cut_min = -20 # cut data before this time 
t_cut_max = 30 # cut data after this time
# Best fit cut-off power law info:
alpha = -1.196 # photon index
ep = 86.386 # keV, peak energy 
norm = 1.32e-01 # counts cm−2 s^−1 keV^−1, normalization