name = "130408A"
z = 3.758 # Measured redshift 
zmax = 10. # Highest redshift to simulate at
fn = "data_files/grb_130408A/grb_130408A_1chan_1s.lc" # file path to light curve
t_true = 4.240 # true T90
t_cut_min = -10 # cut data before this time 
t_cut_max = 20 # cut data after this time
# Best fit cut-off power law info:
alpha = -0.525 # photon index
ep = 147.642 # keV, peak energy 
norm = 6.32e-02 # counts cm−2 s^−1 keV^−1, normalization