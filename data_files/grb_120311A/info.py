name = "120311A"
fn = "data_files/grb_120311A/grb_120311A_1chan_64ms.lc" # file path to light curve

z = 0.350 # Measured redshift 
zmax = 2.5 # Highest redshift to simulate at

t_true = 8.86 # true T90
t_cut_min = -3 # cut data before this time 
t_cut_max = 5 # cut data after this time
t_buffer = 15 # buffer time added before and 
F_true = 0.608 # total mask-weighted counts in T100 

# Best fit cut-off power law info:
alpha = -2.294 # photon index
ep = 150 # keV, peak energy 
norm = 5.82e-02 # counts cm−2 s^−1 keV^−1, normalization