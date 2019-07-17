#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 12:19:06 2019

@author: toma
"""

import numpy as np
import time

#calculates the partition function and write to a file for molecules that has no partition

######################

#store the value of speed of light in cm/s
c = 2.99792458e10
#planck constant erg*s
h = 6.62606885e-27
#boltzmann constant
k_B = 1.38064852e-16
#reference Temperature in K
T_ref = 296 
#store the value of c_2 = h * c / k_B
c_2 = h * c / k_B
#store the conversion from gram to amu
G_TO_AMU = 1.66054e-24 #1.66053904e-24

########################

start_time = time.time()

def calculate_partition(states_filepath, partition_filepath, max_temp): #using HITRAN formula
    Es, gs= np.loadtxt(states_filepath, usecols=(1, 2), unpack=True)
    length_states = len(Es)
    '''
    with open(partition_filepath, 'w') as partitions_file:
        #calculate partiton for temperature from 1 K to ??? K
        for T in range(1, 9001):
            partition_value = 0
            #sum up the partition value for each line
            for i in range(length_states):
                partition_value += gs[i] * np.exp(-c_2 * Es[i] / T)
            partitions_file.write(str(float(T)) + ' ' + str(partition_value) + ' \n')
    partitions_file.close()
    '''
    #significant speed up using numpy array
    #calculate partiton for temperature from 1 K to max_temp K
    Ts = np.arange(1, max_temp + 1)
    partition_values = np.zeros(max_temp)
    for i in range(length_states):
        partition_values += gs[i] * np.exp(-c_2 * Es[i] / Ts)
        #in fact till this step it is good enough...the Ts and partition_values can
        #be passed into exomol_import if we dont care about formating
    output = np.column_stack((Ts.flatten(), partition_values.flatten()))
    #print(output)
    np.savetxt(partition_filepath, output)

########################
    
calculate_partition('/home/toma/Desktop/linelists-database/PH3_states.txt', '/home/toma/Desktop/1314258', 3000)

print("Finished in %s seconds" % (time.time() - start_time))