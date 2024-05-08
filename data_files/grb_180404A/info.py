name = "180404A"
fn = "data_files/grb_180404A/grb_180404A_1chan_1s.lc" # file path to light curve

z = 1.00 # Measured redshift 
zmax =  # Highest redshift to simulate at

t_true = 36.280 # true T90
t_cut_min = # cut data before this time 
t_cut_max = # cut data after this time
t_buffer =  # buffer time added before and 
F_true =  # total mask-weighted counts in T100 

# Best fit cut-off power law info:
alpha = -1.581 # photon index
ep = 202.949 # keV, peak energy 
norm = 1.13e-02 # counts cm−2 s^−1 keV^−1, normalization