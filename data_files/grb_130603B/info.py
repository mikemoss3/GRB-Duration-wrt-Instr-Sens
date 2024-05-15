name = "130603B"
fn = "data_files/grb_130603B/grb_130603B_1chan_64ms.lc" # file path to light curve

z = 0.356 # Measured redshift 
zmax = 1.4 # Highest redshift to simulate at

t_true = 0.19 # true T90
t_start = 0.012 # start true T90
t_cut_min = -2 # cut data before this time 
t_cut_max = 2 # cut data after this time
t_buffer = 20 # buffer time added before and 
F_true = 0.529 # total mask-weighted counts in T100 

# Best fit cut-off power law info:
alpha = -0.750 # photon index
ep = 660 # keV, peak energy 
norm = 1.3e-01 # counts cm−2 s^−1 keV^−1, normalization