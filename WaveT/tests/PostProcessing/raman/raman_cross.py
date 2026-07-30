#!/usr/bin/env python3
import numpy as np
import math

# ATTENTION: change only values in name_list, sigma and field amplitude

name_list = ['raman_time_np_20.dat','raman_time_np_30.dat','raman_time_np_40.dat','raman_time_np_50.dat','raman_time_np_60.dat','raman_time_np_80.dat','raman_time_np_100.dat']

sigma = 10526.82750997          # pulse width in atomic units (sigma in WaveT input)
light_speed = 137.036           # light speed in atomic units
field_amplitude = 10**(-6)      # electric field amplitude in atomic units (fmax in WaveT input)

for i in range(len(name_list)):
    filename=name_list[i]
    raman_read=np.loadtxt(filename)
    r,c=raman_read.shape
    
    cross_section=np.zeros(shape=(r,1))
    raman=np.zeros(shape=(r,2))
    cross_section[:,0]=raman_read[:,1]*2/(math.pi**(3/2)*light_speed**4*sigma*field_amplitude**2)
    
    raman[:,0]=raman_read[:,0]
    raman[:,1]=np.multiply(cross_section[:,0],np.power(raman_read[:,0],4))
    
    np.savetxt(name_list[i], raman, delimiter="\t")   # ATTENTION: new file are overwritten to old ones
