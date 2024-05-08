name = "100816A"
fn = "data_files/grb_100816A/grb_100816A_1chan_64ms.lc" # file path to light curve

z = 0.805 # Measured redshift 
zmax = 5. # Highest redshift to simulate at

t_true = 2.884 # true T90
t_cut_min = -2 # cut data before this time 
t_cut_max = 6 # cut data after this time
t_buffer = 20 # buffer time added before and 
F_true = 2.679 # total mask-weighted counts in T100 

# Best fit cut-off power law info:
alpha = -0.736 # photon index
ep = 172.191 # keV, peak energy 
norm = 1.50e-01 # counts cm−2 s^−1 keV^−1, normalization