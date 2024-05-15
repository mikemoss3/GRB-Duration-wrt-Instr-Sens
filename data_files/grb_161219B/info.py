name = "161219B"
fn = "data_files/grb_161219B/grb_161219B_1chan_64ms.lc" # file path to light curve

z = 0.148 # Measured redshift 
zmax = 0.7 # Highest redshift to simulate at

t_true = 21.80 # true T90
t_start = 0.048 # start true T90
t_cut_min = -2 # cut data before this time 
t_cut_max = 25 # cut data after this time
t_buffer = 30 # buffer time added before and 
F_true = 2.568 # total mask-weighted counts in T100 

# Best fit cut-off power law info:
alpha = -1.294 # photon index
ep = 61.877 # keV, peak energy 
norm = 1.55e-01 # counts cm−2 s^−1 keV^−1, normalization