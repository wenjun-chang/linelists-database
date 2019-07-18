#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 12:27:36 2019

@author: toma
"""

#imports all HITRAN molecule data into the database except for ClONO2, SF6, and CF4.
#only external links to HITRAN 2012 data of those molecules are given
#total molecules = 125 - 3 = 122

import hapi
import numpy as np
import time
import insert_hitran
from query_functions import sql_order
import os

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
mol_ids, iso_ids, iso_names, iso_abundances, iso_masses, mol_names = \
np.loadtxt('/home/toma/Desktop/molecule_properties (copy).txt', dtype='str', skiprows=1, usecols=(1, 2, 3, 4, 5, 6), unpack=True)

for i in range(len(mol_ids)):
    version_name = 'HITRAN_2016'
    
    particle_property_query = "INSERT INTO particles VALUES('%s', '%s', '%s', '%s', '%s', null);" % (mol_names[i], iso_names[i], \
                                                           iso_abundances[i], iso_masses[i], version_name)
    #insert each molecule's properties into particles table
    sql_order(particle_property_query)
    '''
    #then, fetch all the data from HITRAN using HAPI
    hapi.db_begin('data')
    #becasue cannot choose inifinity as upper limit, use a giant number instead
    #gpp is upper state degeneracy
    hapi.fetch(mol_names[i], int(mol_ids[i]), int(iso_ids[i]), 0, 1e9, Parameters=['nu', 'a', 'gamma_air', 'n_air', 'delta_air', \
               'elower', 'gpp', 'gamma_H2', 'n_H2', 'delta_H2', 'gamma_He', 'n_He', 'delta_He'])
    
    #open the file and use insert_hitran.py to insert all parameters into transitions table
    filename = '/home/toma/Desktop/linelists-database/data/{}.data'.format(mol_names[i])
    
    insert_hitran.insert_hitran(filename, version_name, i + 1)
    
    #delete the files since the files are named by HAPI using mol_name instead of iso_name
    #s.t. python wont get confused in the for loop
    header_filename = '/home/toma/Desktop/linelists-database/data/{}.header'.format(mol_names[i])
    os.remove(filename)
    os.remove(header_filename)
    '''
print("Finished in %s seconds" % (time.time() - start_time))