name = "111008A"
z = 4.99 # Measured redshift 
zmax = 10. # Highest redshift to simulate at
fn = "data_files/grb_111008A/grb_111008A_1chan_1s.lc" # file path to light curve
t_true = 62.848 # true T90
F_true = 4.29 # counts / cm^-2 
t_cut_min = -40 # cut data before this time 
t_cut_max = 100 # cut data after this time
# Best fit cut-off power law info:
alpha = -1.698 # photon index
ep = 122.566 # keV, peak energy 
norm = 9.12e-03 # counts cm−2 s^−1 keV^−1, normalization