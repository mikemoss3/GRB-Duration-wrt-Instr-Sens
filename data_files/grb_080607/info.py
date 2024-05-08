name = "080607"
fn = "data_files/grb_080607/grb_080607_1chan_1s.lc" # file path to light curve

z = 3.04 # Measured redshift 


t_true = 78.972 # true T90
t_cut_min = -7 # cut data before this time 
t_cut_max = 160 # cut data after this time
t_buffer = 70 # buffer time added before and 
F_true = 35.144 # total mask-weighted counts in T100 

# Best fit cut-off power law info:
alpha = -0.659 # photon index
ep = 552.387 # keV, peak energy 
norm = 1.35e-01 # counts cm−2 s^−1 keV^−1, normalization