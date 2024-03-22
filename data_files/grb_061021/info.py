name = "061021"
z = 0.35 # Measured redshift 
zmax = 4. # Highest redshift to simulate at
fn = "data_files/grb_061021/grb_061021_1chan_1s.lc" # file path to light curve
t_true = 47.820 # true T90
t_cut_min = -60 # cut data before this time 
t_cut_max = 100 # cut data after this time
# Best fit cut-off power law info:
alpha = -1.269 # photon index
ep = 9973.470 # keV, peak energy 
norm = 9.74e-02 # counts cm−2 s^−1 keV^−1, normalization