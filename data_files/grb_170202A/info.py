name = "170202A"
z = 3.645 # Measured redshift 
zmax = 1. # Highest redshift to simulate at
fn = "data_files/grb_170202A/grb_170202A_1chan_1s.lc" # file path to light curve
t_true = 37.760 # true T90
t_cut_min = -30 # cut data before this time 
t_cut_max = 100 # cut data after this time
# Best fit cut-off power law info:
alpha = -1.413 # photon index
ep = 113.258 # keV, peak energy 
norm = 6.35e-03 # counts cm−2 s^−1 keV^−1, normalization