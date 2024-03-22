name = "100816A"
z = 0.805 # Measured redshift 
zmax = 5. # Highest redshift to simulate at
fn = "data_files/grb_100816A/grb_100816A_1chan_64ms.lc" # file path to light curve
t_true = 2.884 # true T90
t_cut_min = -10 # cut data before this time 
t_cut_max = 20 # cut data after this time
# Best fit cut-off power law info:
alpha = -0.736 # photon index
ep = 172.191 # keV, peak energy 
norm = 1.50e-01 # counts cm−2 s^−1 keV^−1, normalization