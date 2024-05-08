name = "060206"
fn = "data_files/grb_060206/grb_060206_1chan_64ms.lc" # file path to light curve

z = 4.05 # Measured redshift 


t_true = 7.552 # true T90
t_cut_min = -2 # cut data before this time 
t_cut_max = 15 # cut data after this time
t_buffer = 22 # buffer time added before and 
F_true = 1.553 # total mask-weighted counts in T100 

# Best fit cut-off power law info:
alpha = -0.909 # photon index
ep = 81.764 # keV, peak energy 
norm = 1.67e-01 # counts cm−2 s^−1 keV^−1, normalization