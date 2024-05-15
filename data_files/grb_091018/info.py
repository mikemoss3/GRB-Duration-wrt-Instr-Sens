name = "091018"
fn = "data_files/grb_091018/grb_091018_1chan_64ms.lc" # file path to light curve

z = 0.97 # Measured redshift 
zmax = 4.5 # Highest redshift to simulate at

t_true = 6.368 # true T90
t_start = 0.240 # start true T90
t_cut_min = -2 # cut data before this time 
t_cut_max = 10 # cut data after this time
t_buffer = 30 # buffer time added before and 
F_true = 3.086 # total mask-weighted counts in T100 

# Best fit cut-off power law info:
alpha = -1.758 # photon index
ep = 19.425 # keV, peak energy 
norm = 2.20e-01 # counts cm−2 s^−1 keV^−1, normalization