name = "091127"
fn = "data_files/grb_091127/grb_091127_1chan_64ms.lc" # file path to light curve

z = 0.49 # Measured redshift 
zmax = 5. # Highest redshift to simulate at

t_true = 9.35 # true T90
t_start = -0.116 # start true T90
t_cut_min = -2 # cut data before this time 
t_cut_max = 12 # cut data after this time
t_buffer = 30 # buffer time added before and 
F_true = 13.086 # total mask-weighted counts in T100 

# Best fit cut-off power law info:
alpha = -1.797 # photon index
ep = 46.374 # keV, peak energy 
norm = 1.18e-01 # counts cm−2 s^−1 keV^−1, normalization