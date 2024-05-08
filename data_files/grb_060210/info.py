name = "060210"
fn = "data_files/grb_060210/grb_060210_1chan_1s.lc" # file path to light curve

z = 3.91 # Measured redshift 


t_true = 288.000 # true T90
t_cut_min = -240 # cut data before this time 
t_cut_max = 140 # cut data after this time
t_buffer = 100 # buffer time added before and 
F_true = 11.4 # total mask-weighted counts in T100 

# Best fit cut-off power law info:
alpha = -1.186 # photon index
ep = 150.080  # keV, peak energy 
norm = 1.22e-01 # counts cm−2 s^−1 keV^−1, normalization