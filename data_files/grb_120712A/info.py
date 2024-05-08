name = "120712A"
fn = "data_files/grb_120712A/grb_120712A_1chan_1s.lc" # file path to light curve

z = 4.175 # Measured redshift 


t_true = 14.808 # true T90
t_cut_min = -6 # cut data before this time 
t_cut_max = 18 # cut data after this time
t_buffer = 20 # buffer time added before and 
F_true = 2.680 # total mask-weighted counts in T100 

# Best fit cut-off power law info:
alpha = -0.989 # photon index
ep = 144.267 # keV, peak energy 
norm = 1.35e-02 # counts cm−2 s^−1 keV^−1, normalization