name = "160804A"
fn = "data_files/grb_160804A/grb_160804A_1chan_1s.lc" # file path to light curve

z = 0.736 # Measured redshift 
zmax = 4. # Highest redshift to simulate at

t_true = 185.87 # true T90
t_start = 38.204 # start true T90
t_cut_min = -50 # cut data before this time 
t_cut_max = 300 # cut data after this time
t_buffer = 200 # buffer time added before and 
F_true = 19.294 # total mask-weighted counts in T100 

# Best fit cut-off power law info:
alpha = -1.415 # photon index
ep = 78.961 # keV, peak energy 
norm = 1.3e-01 # counts cm−2 s^−1 keV^−1, normalization