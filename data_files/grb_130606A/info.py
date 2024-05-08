name = "130606A"
fn = "data_files/grb_130606A/grb_130606A_1chan_1s.lc" # file path to light curve
z = 5.91 # Measured redshift 


t_true = 276.656 # true T90
t_cut_min = -2 # cut data before this time 
t_cut_max = 300 # cut data after this time
t_buffer = 100 # buffer time added before and 
F_true = 4.701 # total mask-weighted counts in T100 

# Best fit cut-off power law info:
alpha = -1.291 # photon index
ep = 150.299 # keV, peak energy 
norm = 1.28e-03 # counts cm−2 s^−1 keV^−1, normalization