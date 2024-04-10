name = "120712A"
z = 4.175 # Measured redshift 
zmax = 10. # Highest redshift to simulate at
fn = "data_files/grb_120712A/grb_120712A_1chan_1s.lc" # file path to light curve
t_true = 14.808 # true T90
F_true = 2.684 # counts / cm^-2 
t_cut_min = -20 # cut data before this time 
t_cut_max = 40 # cut data after this time
# Best fit cut-off power law info:
alpha = -0.989 # photon index
ep = 144.267 # keV, peak energy 
norm = 1.35e-02 # counts cm−2 s^−1 keV^−1, normalization