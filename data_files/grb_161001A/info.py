name = "161001A"
fn = "data_files/grb_161001A/grb_161001A_1chan_64ms.lc" # file path to light curve

z = 0.891 # Measured redshift 
zmax = 4.5 # Highest redshift to simulate at

t_true = 12.20 # true T90
t_cut_min = -2 # cut data before this time 
t_cut_max = 5 # cut data after this time
t_buffer = 15 # buffer time added before and 
F_true = 0.955 # total mask-weighted counts in T100 

# Best fit cut-off power law info:
alpha = -0.841 # photon index
ep = 282.469 # keV, peak energy 
norm = 1.34e-01 # counts cm−2 s^−1 keV^−1, normalization