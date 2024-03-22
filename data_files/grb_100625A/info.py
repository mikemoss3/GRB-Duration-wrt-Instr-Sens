name = "100625A"
z = 0.452 # Measured redshift 
zmax = 2.0 # Highest redshift to simulate at
fn = "data_files/grb_100625A/grb_100625A_1chan_64ms.lc" # file path to light curve
t_true = 0.332 # true T90
t_cut_min = -10 # cut data before this time 
t_cut_max = 20 # cut data after this time
# Best fit cut-off power law info:
alpha = -0.724 # photon index
ep = 446.613 # keV, peak energy 
norm = 1.35e-01 # counts cm−2 s^−1 keV^−1, normalization