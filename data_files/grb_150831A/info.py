name = "150831A"
fn = "data_files/grb_150831A/grb_150831A_1chan_64ms.lc" # file path to light curve

z = 1.18 # Measured redshift 
zmax =  # Highest redshift to simulate at

t_true = 0.920 # true T90
t_cut_min =  # cut data before this time 
t_cut_max = # cut data after this time
t_buffer =  # buffer time added before and 
F_true = # total mask-weighted counts in T100 

# Best fit cut-off power law info:
alpha = -0.712 # photon index
ep = 9999.360 # keV, peak energy 
norm = 2.06e-02  # counts cm−2 s^−1 keV^−1, normalization