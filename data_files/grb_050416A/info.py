name = "050416A"
z = 0.65 # Measured redshift 
zmax = 2. # Highest redshift to simulate at
fn = "data_files/grb_050416A/grb_050416A_1chan_64ms.lc" # file path to light curve
t_true = 4.2 # true T90
t_cut_min = -20 # cut data before this time 
t_cut_max = 30 # cut data after this time
# Best fit cut-off power law info:
alpha = 0.097 # photon index
ep = 22.692 # keV, peak energy 
norm = 8.46e-1 # counts cm−2 s^−1 keV^−1, normalization