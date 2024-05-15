name = "110715A"
fn = "data_files/grb_110715A/grb_110715A_1chan_1s.lc" # file path to light curve

z = 0.823 # Measured redshift 
zmax = 6. # Highest redshift to simulate at

t_true = 13.0 # true T90
t_start = -0.144 # start true T90
t_cut_min = -5 # cut data before this time 
t_cut_max = 20 # cut data after this time
t_buffer = 50 # buffer time added before and 
F_true = 18.7 # total mask-weighted counts in T100 

# Best fit cut-off power law info:
alpha = -1.254 # photon index
ep = 119.804 # keV, peak energy 
norm = 1.3e-01 # counts cm−2 s^−1 keV^−1, normalization