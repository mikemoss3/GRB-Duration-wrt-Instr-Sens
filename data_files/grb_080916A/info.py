name = "080916A"
fn = "data_files/grb_080916A/grb_080916A_1chan_1s.lc" # file path to light curve

z = 0.69 # Measured redshift 
zmax = 3. # Highest redshift to simulate at

t_true = 61.348 # true T90
t_cut_min = -4 # cut data before this time 
t_cut_max = 90 # cut data after this time
t_buffer = 50 # buffer time added before and 
F_true = 7.041 # total mask-weighted counts in T100 

# Best fit cut-off power law info:
alpha = -0.164 # photon index
ep = 138.254 # keV, peak energy 
norm = 2.02e-01 # counts cm−2 s^−1 keV^−1, normalization