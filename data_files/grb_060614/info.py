name = "060614"
fn = "data_files/grb_060614/grb_060614_1chan_1s.lc" # file path to light curve

z = 0.1254 # Measured redshift 
zmax = 1. # Highest redshift to simulate at

t_true = 109.104 # true T90
t_start = 0.332 # start true T90
t_cut_min = -3 # cut data before this time 
t_cut_max = 120 # cut data after this time
t_buffer = 200 # buffer time added before and 
F_true = 38.715 # total mask-weighted counts in T100 

# Best fit cut-off power law info:
alpha = -1.382 # photon index
ep = 147.910 # keV, peak energy 
norm = 1.15e-01 # counts cm−2 s^−1 keV^−1, normalization