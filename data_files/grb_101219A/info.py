name = "101219A"
fn = "data_files/grb_101219A/grb_101219A_1chan_64ms.lc" # file path to light curve

z = 0.718 # Measured redshift 
zmax = 2 # Highest redshift to simulate at

t_true = 2.17 # true T90
t_start = 0.004 # start true T90
t_cut_min = -2 # cut data before this time 
t_cut_max = 4 # cut data after this time
t_buffer = 20 # buffer time added before and 
F_true = 0.534 # total mask-weighted counts in T100 

# Best fit cut-off power law info:
alpha = -0.621 # photon index
ep = 490 # keV, peak energy 
norm = 1.4e-01 # counts cm−2 s^−1 keV^−1, normalization