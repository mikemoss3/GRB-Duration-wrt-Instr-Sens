name = "100625A"
fn = "data_files/grb_100625A/grb_100625A_1chan_64ms.lc" # file path to light curve

z = 0.452 # Measured redshift 
zmax = 2.0 # Highest redshift to simulate at

t_true = 16.34 # true T90
t_cut_min = -2 # cut data before this time 
t_cut_max = 2 # cut data after this time
t_buffer = 12 # buffer time added before and 
F_true = 0.284 # total mask-weighted counts in T100 

# Best fit cut-off power law info:
alpha = -0.724 # photon index
ep = 446.613 # keV, peak energy 
norm = 1.35e-01 # counts cm−2 s^−1 keV^−1, normalization