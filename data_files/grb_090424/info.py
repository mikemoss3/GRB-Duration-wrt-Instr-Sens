name = "090424"
fn = "data_files/grb_090424/grb_090424_1chan_1s.lc" # file path to light curve

z = 0.54 # Measured redshift 
zmax = 7 # Highest redshift to simulate at

t_true = 49.460 # true T90
t_cut_min = -2 # cut data before this time 
t_cut_max = 103 # cut data after this time
t_buffer = 70 # buffer time added before and 
F_true = 33. # total mask-weighted counts in T100 

# Best fit cut-off power law info:
alpha = -0.984 # photon index
ep = 183.397 # keV, peak energy 
norm = 1.26e-01 # counts cm−2 s^−1 keV^−1, normalization