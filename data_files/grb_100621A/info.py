name = "100621A"
z = 0.54 # Measured redshift 
zmax = 8.5 # Highest redshift to simulate at
fn = "data_files/grb_100621A/grb_100621A_1chan_64ms.lc" # file path to light curve
t_true = 63.552 # true T90
t_cut_min = -70 # cut data before this time 
t_cut_max = 260 # cut data after this time
# Best fit cut-off power law info:
alpha = -1.812 # photon index
ep = 129.012 # keV, peak energy 
norm = 8.55e-02 # counts cm−2 s^−1 keV^−1, normalization