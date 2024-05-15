name = "080430"
fn = "data_files/grb_080430/grb_080430_1chan_1s.lc" # file path to light curve

z = 0.77 # Measured redshift 
zmax = 2 # Highest redshift to simulate at

t_true = 13.872 # true T90
t_start = -0.3 # start true T90
t_cut_min = -3 # cut data before this time 
t_cut_max = 20 # cut data after this time
t_buffer = 100 # buffer time added before and after emission
F_true = 2.146 # total mask-weighted counts in T100 

# Best fit cut-off power law info:
alpha = -1.645 # photon index
ep = 151.58 # keV, peak energy 
norm = 2.25e-02 # counts cm−2 s^−1 keV^−1, normalization
