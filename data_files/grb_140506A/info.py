name = "140506A"
fn = "data_files/grb_140506A/grb_140506A_1chan_1s.lc" # file path to light curve

z = 0.889 # Measured redshift 
zmax = 7. # Highest redshift to simulate at

t_true = 111.104 # true T90
t_cut_min = -2 # cut data before this time 
t_cut_max = 125 # cut data after this time
t_buffer = 60 # buffer time added before and 
F_true = 4.079 # total mask-weighted counts in T100 

# Best fit cut-off power law info:
alpha = -1.449 # photon index
ep = 271.007 # keV, peak energy 
norm = 9.65e-02 # counts cm−2 s^−1 keV^−1, normalization