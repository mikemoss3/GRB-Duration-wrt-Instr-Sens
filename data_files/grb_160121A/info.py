name = "160121A"
fn = "data_files/grb_160121A/grb_160121A_1chan_64ms.lc" # file path to light curve

z = 1.960 # Measured redshift 
zmax =  # Highest redshift to simulate at

t_true = 10.500 # true T90
t_cut_min =  # cut data before this time 
t_cut_max = # cut data after this time
t_buffer =  # buffer time added before and 
F_true =  # total mask-weighted counts in T100 

# Best fit cut-off power law info:
alpha = -0.929 # photon index
ep = 51.480  # keV, peak energy 
norm = 2.74e-02 # counts cm−2 s^−1 keV^−1, normalization