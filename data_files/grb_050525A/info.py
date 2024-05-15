name = "050525A"
fn = "data_files/grb_050525A/grb_050525A_1chan_64ms.lc" # file path to light curve

z = 0.61 # Measured redshift 
zmax = 6.5 # Highest redshift to simulate at

t_true = 8.836 # true T90
t_start = 0.328 # start true T90
t_cut_min = -2 # cut data before this time 
t_cut_max = 14 # cut data after this time
t_buffer = 50 # buffer time added before and 
F_true = 25.613  # total mask-weighted counts in T100 

# Best fit cut-off power law info:
alpha = -1.026 # photon index
ep = 80.4 # keV, peak energy 
norm = 1.67e-01  # counts cm−2 s^−1 keV^−1, normalization