name = "111228A"
fn = "data_files/grb_111228A/grb_111228A_1chan_1s.lc" # file path to light curve

z = 0.716 # Measured redshift 
zmax = 7. # Highest redshift to simulate at

t_true = 101.244 # true T90
t_cut_min = -15 # cut data before this time 
t_cut_max = 120 # cut data after this time
t_buffer = 130 # buffer time added before and 
F_true = 16.025 # total mask-weighted counts in T100 

# Best fit cut-off power law info:
alpha = -1.644 # photon index
ep = 93.714 # keV, peak energy 
norm = 1.07e-01 # counts cm−2 s^−1 keV^−1, normalization