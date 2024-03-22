name = "161001A"
z = 0.891 # Measured redshift 
zmax = 4. # Highest redshift to simulate at
fn = "data_files/grb_161001A/grb_161001A_1chan_64ms.lc" # file path to light curve
t_true = 2.6 # true T90
t_cut_min = -10 # cut data before this time 
t_cut_max = 20 # cut data after this time
# Best fit cut-off power law info:
alpha = -0.841 # photon index
ep = 282.469 # keV, peak energy 
norm = 1.34e-01 # counts cm−2 s^−1 keV^−1, normalization