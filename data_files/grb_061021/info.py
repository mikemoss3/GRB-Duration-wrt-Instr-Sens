name = "061021"
fn = "data_files/grb_061021/grb_061021_1chan_1s.lc" # file path to light curve

z = 0.35 # Measured redshift 
zmax = 4. # Highest redshift to simulate at

t_true = 56.12 # true T90
t_cut_min = -2 # cut data before this time 
t_cut_max = 66 # cut data after this time
t_buffer = 40 # buffer time added before and 
F_true = 4.480 # total mask-weighted counts in T100 

# Best fit cut-off power law info:
alpha = -1.269 # photon index
ep = 777 # keV, peak energy 
norm = 1.02e-01 # counts cm−2 s^−1 keV^−1, normalization