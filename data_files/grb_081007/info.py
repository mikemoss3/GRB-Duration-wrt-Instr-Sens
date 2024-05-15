name = "081007"
fn = "data_files/grb_081007/grb_081007_1chan_64ms.lc" # file path to light curve

z = 0.5295 # Measured redshift 
zmax = 2. # Highest redshift to simulate at

t_true = 12.5 # true T90
t_start = -4 # start true T90
t_cut_min = -8 # cut data before this time 
t_cut_max = 20 # cut data after this time
t_buffer = 20 # buffer time added before and 
F_true = 1.336 # total mask-weighted counts in T100 

# Best fit cut-off power law info:
alpha = -1.52043 # photon index
ep = 20.1229 # keV, peak energy 
norm = 2.77485E-01 # counts cm−2 s^−1 keV^−1, normalization