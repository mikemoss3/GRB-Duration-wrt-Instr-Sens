name = "170202A"
fn = "data_files/grb_170202A/grb_170202A_1chan_1s.lc" # file path to light curve

z = 3.645 # Measured redshift 

t_true = 37.760 # true T90
t_cut_min = -1 # cut data before this time 
t_cut_max = 80 # cut data after this time
t_buffer = 100 # buffer time added before and 
F_true = 5.595 # total mask-weighted counts in T100 

# Best fit cut-off power law info:
alpha = -1.413 # photon index
ep = 113.258 # keV, peak energy 
norm = 6.35e-03 # counts cm−2 s^−1 keV^−1, normalization