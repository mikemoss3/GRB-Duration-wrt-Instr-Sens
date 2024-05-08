name = "060927"
fn = "data_files/grb_060927/grb_060927_1chan_1s.lc" # file path to light curve

z = 5.47 # Measured redshift 


t_true = 22.416 # true T90
t_cut_min = -2 # cut data before this time 
t_cut_max = 25 # cut data after this time
t_buffer = 30 # buffer time added before and 
F_true = 2. # total mask-weighted counts in T100 

# Best fit cut-off power law info:
alpha = -0.808 # photon index
ep = 70.673 # keV, peak energy 
norm = 1.96e-01 # counts cm−2 s^−1 keV^−1, normalization