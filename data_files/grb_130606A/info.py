name = "130606A"
z = 5.91 # Measured redshift 
zmax = 10. # Highest redshift to simulate at
fn = "data_files/grb_130606A/grb_130606A_1chan_1s.lc" # file path to light curve
t_true = 276.656 # true T90
F_true = 4.7 # counts / cm^-2 
t_cut_min = -70 # cut data before this time 
t_cut_max = 350 # cut data after this time
# Best fit cut-off power law info:
alpha = -1.291 # photon index
ep = 150.299 # keV, peak energy 
norm = 1.28e-03 # counts cm−2 s^−1 keV^−1, normalization