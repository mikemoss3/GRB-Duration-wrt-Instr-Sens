name = "060306"
fn = "data_files/grb_060306/grb_060306_1chan_1s.lc" # file path to light curve

z = 3.5 # Measured redshift 


t_true = 60.940 # true T90
t_cut_min = -3 # cut data before this time 
t_cut_max = 70 # cut data after this time
t_buffer = 40 # buffer time added before and 
F_true = 3.549 # total mask-weighted counts in T100 

# Best fit cut-off power law info:
alpha = -1.254 # photon index
ep = 69.380 # keV, peak energy 
norm = 1.50e-01 # counts cm−2 s^−1 keV^−1, normalization