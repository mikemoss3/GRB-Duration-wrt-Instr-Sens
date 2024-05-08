name = "160425A"
fn = "data_files/grb_160425A/grb_160425A_1chan_1s.lc" # file path to light curve

z = 0.555 # Measured redshift 
zmax = 2.6 # Highest redshift to simulate at

t_true = 125.69 # true T90
t_cut_min = -2 # cut data before this time 
t_cut_max = 123 # cut data after this time
t_buffer = 80 # buffer time added before and 
F_true = 4.079 # total mask-weighted counts in T100 

# Best fit cut-off power law info:
alpha = -1.258 # photon index
ep = 129.068 # keV, peak energy 
norm = 1.24e-01 # counts cm−2 s^−1 keV^−1, normalization