name = "060912A"
fn = "data_files/grb_060912A/grb_060912A_1chan_64ms.lc" # file path to light curve

z = 0.94 # Measured redshift 
zmax = 3 # Highest redshift to simulate at

t_true = 11.5 # true T90
t_start = -0.116 # start true T90
t_cut_min = -2 # cut data before this time 
t_cut_max = 17 # cut data after this time
t_buffer = 40 # buffer time added before and 
F_true = 2.25 # total mask-weighted counts in T100 

# Best fit cut-off power law info:
alpha = -1.196 # photon index
ep = 86.386 # keV, peak energy 
norm = 1.32e-01 # counts cm−2 s^−1 keV^−1, normalization