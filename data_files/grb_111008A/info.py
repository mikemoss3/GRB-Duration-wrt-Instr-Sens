name = "111008A"
fn = "data_files/grb_111008A/grb_111008A_1chan_1s.lc" # file path to light curve

z = 4.99 # Measured redshift 


t_true = 62.848 # true T90
t_cut_min = -4 # cut data before this time 
t_cut_max = 70 # cut data after this time
t_buffer = 80 # buffer time added before and 
F_true = 8.492 # total mask-weighted counts in T100 

# Best fit cut-off power law info:
alpha = -1.698 # photon index
ep = 122.566 # keV, peak energy 
norm = 9.12e-03 # counts cm−2 s^−1 keV^−1, normalization