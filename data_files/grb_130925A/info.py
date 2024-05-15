name = "130925A"
fn = "data_files/grb_130925A/grb_130925A_1chan_1s.lc" # file path to light curve

z = 0.347 # Measured redshift 
zmax = 2. # Highest redshift to simulate at

t_true = 285.73 # true T90
t_start = -48.892 # start true T90
t_cut_min = -62 # cut data before this time 
t_cut_max = 400 # cut data after this time
t_buffer = 200 # buffer time added before and 
F_true = 67 # total mask-weighted counts in T100 

# Best fit cut-off power law info:
alpha = -1.847 # photon index
ep = 39.703 # keV, peak energy 
norm = 1.23e-01 # counts cm−2 s^−1 keV^−1, normalization