name = "131004A"
fn = "data_files/grb_131004A/grb_131004A_1chan_64ms.lc" # file path to light curve

z = 0.717 # Measured redshift 
zmax = # Highest redshift to simulate at

t_true = 1.536 # true T90
t_cut_min = -50 # cut data before this time 
t_cut_max = 50 # cut data after this time
t_buffer = 200 # buffer time added before and 
F_true = 520.680 # total mask-weighted counts in T100 

# Best fit cut-off power law info:
alpha = -0.852 # photon index
ep = 63.374 # keV, peak energy 
norm = 3.38e-02 # counts cm−2 s^−1 keV^−1, normalization