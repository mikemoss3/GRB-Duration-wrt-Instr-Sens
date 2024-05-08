name = "050416A"
fn = "data_files/grb_050416A/grb_050416A_1chan_64ms.lc" # file path to light curve

z = 0.65 # Measured redshift 
zmax = 2. # Highest redshift to simulate at

t_true = 4.2 # true T90
t_cut_min = -2 # cut data before this time 
t_cut_max = 11 # cut data after this time
t_buffer = 20 # buffer time added before and 
F_true = 0.925 # total mask-weighted counts in T100 

# Best fit cut-off power law info:
alpha = 0.097 # photon index
ep = 22.692 # keV, peak energy 
norm = 8.46e-1 # counts cm−2 s^−1 keV^−1, normalization