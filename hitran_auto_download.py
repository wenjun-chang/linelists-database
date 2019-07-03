#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 12:27:36 2019

@author: toma
"""

import hapi
import numpy as np
import time
import insert_hitran

#this thing gives us the table of molecule properties in shell
#saved it manually in a file called /home/toma/Desktop/molecule_properties.txt
#getHelp(ISO_ID)

start_time = time.time()

'''
#change file to loadable version
with open('/home/toma/Desktop/molecule_properties.txt', 'r') as f: 
    outfile = open('/home/toma/Desktop/molecule_properties (copy).txt', 'w')
    for line in f: 
        data = line.split()
        data.pop(1)
        for i in data:
            outfile.write("%s " % i)
        outfile.write('\n')
    outfile.close()
f.close()
'''    

#everything is in string though to be noticed
#length of lists are 124
mol_id, iso_id, iso_name, iso_abundance, iso_mass, mol_name = \
np.loadtxt('/home/toma/Desktop/molecule_properties (copy).txt', dtype='str', skiprows=1, usecols=(1, 2, 3, 4, 5, 6), unpack=True)

#don't know why this not workkkkkkkkkkkkkkkkkkkk???????????
hapi.db_begin('data')
hapi.fetch('CO', 5, 1, 0, 100000000000000000, ParameterGroups=['Standard', 'Voigt_Air', 'Voigt_H2', 'Voigt_He'], Parameters=['nu',\
           'a', 'gamma_air', 'n_air', 'delta_air', 'elower', 'gp', 'gamma_H2', 'n_H2', 'delta_H2', 'gamma_He', 'n_He', 'delta_He'])
    
#if above not work for formatting use this to select params from file
#hapi.select('CO', ParameterNames=('nu', 'a', 'gamma_air', 'n_air', 'delta_air', 'elower', 'gp', 'gamma_H2', 'n_H2', 'delta_H2', 'gamma_He', 'n_He', 'delta_He'))

    
filename = '/home/toma/Desktop/linelists-database/data/##.##'
#ideally then know the filename and import the data
#insert_hitran.insert_hitran(filename)
print("Finished in %s seconds" % (time.time() - start_time))