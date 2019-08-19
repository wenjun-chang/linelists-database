#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 17:47:11 2019

@author: toma
"""

#this file get the partition files for hitran on its website and insert into the database

import time
import numpy as np
from query_functions import fetch
from exomol_auto_download import download_file
from exomol_import import insert_partitions

def insert_partition_for_hitran():
    
    start_time = time.time() 
    global_ids, iso_names = np.loadtxt('/home/toma/Desktop/molecule_properties (copy).txt', dtype=str, skiprows=1, usecols=(0, 3), unpack=True)    
    for i in range(len(global_ids)):

        global_id = global_ids[i]
        iso_name = iso_names[i]
        
        if iso_name == '(16O)':
            continue #since it has no partitons file on HITRAN website     
            
        #get particle_id
        get_particle_id = "SELECT particle_id FROM particles WHERE iso_name = '{}'".format(iso_name)
        data = fetch(get_particle_id)
        if data == (): 
            print('empty moelcule encountered')
            continue #note for ClONO2, SF6 and CF4 #because only external links to HITRAN 2012 data of those molecules are given
        particle_id = data[0][0]
        
        #get hitran line source id for that isotopologue
        get_line_sources = "SELECT line_source, line_source_id FROM source_properties WHERE particle_id = {}".format(particle_id)
        sources = fetch(get_line_sources)
        hitran_id = -1
        for one_source in sources: 
            source_name = one_source[0]
            if source_name.startswith('HITRAN'):
                hitran_id = one_source[1]
        if hitran_id == -1:
            raise Exception('Oh Damn this isotopologue is in HITRAN but has no HITRAN lines in the database~')
        
        print(global_id, iso_name, particle_id, hitran_id)
        #download partition file
        
        url = 'https://hitran.org/data/Q/q{}.txt'.format(global_id)
        filename = iso_name + '_hitran_partiitons'
        
        download_file(url, filename)
        download_file(url, filename + '_modified')
        
        #modify the file becasue some of the files have the same problem as the hitemp ones, merging two columns
        modified_file = open('/home/toma/Desktop/linelists-database/' + filename + '_modified', 'w')
        with open('/home/toma/Desktop/linelists-database/' + filename) as file:
            for line in file: 
                T = line[:4]
                par = line[4:]
                modified_file.write(T + ' ' + par + '\n')
        modified_file.close()
        file.close()
        
        #insert the partition file
        insert_partitions('/home/toma/Desktop/linelists-database/' + filename + '_modified', hitran_id, particle_id)
        
        print('Finished inserting hitran partitions for ' + iso_name)
    
    print("Finished in %s seconds" % (time.time() - start_time))
    
if __name__ == '__main__':
    insert_partition_for_hitran()
        
        
        
