name = "150314A"
fn = "data_files/grb_150314A/grb_150314A_1chan_64ms.lc" # file path to light curve

z = 1.758 # Measured redshift 
zmax =  # Highest redshift to simulate at

t_true = 14.780 # true T90
t_cut_min =  # cut data before this time 
t_cut_max =  # cut data after this time
t_buffer = 200 # buffer time added before and 
F_true = # total mask-weighted counts in T100 

# Best fit cut-off power law info:
alpha = -0.465417 # photon index
ep = 248.872 # keV, peak energy 
norm = 0.0262085  # counts cm−2 s^−1 keV^−1, normalization