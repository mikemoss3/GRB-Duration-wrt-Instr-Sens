name = "130427A_cut"
fn = "data_files/grb_130427A_cut/grb_130427A_cut_1chan_1s.lc" # file path to light curve

z = 0.340 # Measured redshift 
zmax = 10. # Highest redshift to simulate at

t_true = 50 # true T90
t_cut_min = -60 # cut data before this time 
t_cut_max = 0 # cut data after this time
t_buffer = 50 # buffer time added before and 
F_true = 392.979 # total mask-weighted counts in T100 

# Best fit cut-off power law info:
alpha = -1.169 # photon index
ep = 824.9 # keV, peak energy 
norm = 1.07e-01 # counts cm−2 s^−1 keV^−1, normalization


