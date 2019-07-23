#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 14:11:38 2019

@author: toma
"""

#this file imports hitemp data

import numpy as np
import MySQLdb
from query_functions import sql_order, fetch
import time
from compute_absorption import get_particle
import zipfile
import os

#############################

#hitemp page: https://hitran.org/hitemp/#ref1
#so far only H2O, CO2, N2O, CO, NO, NO2, OH

#############################

def write_and_sort_hitemp_data(filename, M_id, isotop_num):
    filepath_without_isotop_name = '/home/toma/Desktop/linelists-database/hitemp_'
    print(M_id, isotop_num)
    #get iso name for each iso num
    iso_filepaths = []
    isotop_names = []
    for i in range(1, isotop_num + 1):
        isotop_name = ''
        for j in range(len(mol_ids)):
            if mol_ids[j] == str(M_id) and iso_ids[j] == str(i):
                print('oof')
                isotop_name = iso_names[j]
        isotop_names.append(isotop_name)
        iso_filepaths.append(filepath_without_isotop_name + isotop_name)
    print(isotop_names)
    #read through all lines     
    with open(filename, 'r') as hitemp: 
        for line in hitemp: 
            iso_num = int(line[2:3])
            iso_fp = iso_filepaths[iso_num - 1]
            
            nu = line[3:15].strip()
            A = line[25:35].strip()
            gamma_air = line[35:40].strip()
            elower = line[45:55].strip()
            n_air = line[55:59].strip()
            delta_air = line[59:67].strip()
            g_upper = line[-7:].strip()
            
            #line table arrangement corresponding to tuple indexes: 
            #(nu, A, gamma_air, n_air, delta_air, elower, g_upper, gamma_H2, n_H2, delta_H2, gamma_He, n_He, delta_He, line_source_id, particle_id, line_id)
            #( 0, 1,      2,       3,       4,        5,    6,    7,       8,      9,       10,      11,    12,        13,         14,         15  )
            
            data = [nu, A, gamma_air, n_air, delta_air, elower, g_upper]
            for i in range(len(data)): 
                if data[i].startswith('.'):
                    data[i] = '0' + data[i]
                elif data[i].startswith('-.'):
                    data[i] = '-0' + data[i][1:]
                    
            with open(iso_fp, 'a') as iso_file:  
                for element in data: 
                    iso_file.write("%s " % element)
                iso_file.write("\n")
            iso_file.close()
    hitemp.close()
    return iso_filepaths, isotop_names

##########################
   
def insert_hitemp(fp, isotop_name, line_source_id, ref_link): 
    
    insert_time = time.time()
    
    #connect to the database
    db = MySQLdb.connect(host='localhost', user='toma', passwd='Happy810@', db='linelist')
    
    #create a cursor object
    cursor = db.cursor()
    
    #disable autocommit to improve performance
    sql_order('SET autocommit = 0')
    sql_order('SET unique_checks = 0')
    sql_order('SET foreign_key_checks = 0')
    sql_order('SET sql_log_bin = 0')
    
    try:
        file_length = sum(1 for line in open(fp))
        
        #get particle id
        particle_id = get_particle(isotop_name)[0]
        
        print("Bulk inserting hitemp data...")
        cursor.execute("LOAD DATA LOCAL INFILE '{}' INTO TABLE transitions FIELDS TERMINATED BY ' ' LINES TERMINATED BY '\n' \
                       (@col1, @col2, @col3, @col4, @col5, @col6, @col7) SET nu=@col1, A=@col2, gamma_air=@col3, n_air=@col4, \
                       delta_air=@col5, elower=@col6, g_upper=@col7, line_source_id={}, particle_id={};".format(fp, \
                       line_source_id, particle_id))
        
        #commit changes
        db.commit()
        
        #turn it back on
        sql_order('SET unique_checks = 1')
        sql_order('SET foreign_key_checks = 1')
        sql_order('SET sql_log_bin = 1')
        
        print('Executed {} lines of hitemp data'.format(file_length))
        
    except Exception as e:
        #if errors occur
        db.rollback()
        print('insert hitemp data failed', e)
      
    finally: 
        #close up cursor and connection
        cursor.close()
        db.close()
        
    print("Finished in %s seconds" % (time.time() - insert_time))
    
#########################
    
def insert_all(fp, filename):
    t = time.time()
    iso_num_dict = {'01' : 6, '02' : 7, '04' : 5, '05' : 6, '08' : 3, '10' : 1, '13' : 3}
    M_id = filename[:2]
    isotop_num = iso_num_dict.get(M_id)
    
    #pretty dumb way of splitting verions
    if M_id == '01' or M_id == '02' or M_id == '13':
        line_source_id = hitemp_ids[0]
    else: 
        line_source_id = hitemp_ids[1]
    
    if M_id.startswith('0'):
        M_id = M_id[1:]
        
    iso_filepaths, isotop_names = write_and_sort_hitemp_data(fp, M_id, isotop_num)
    for i in range(len(iso_filepaths)):
        if os.path.isfile(iso_filepaths[i]) is True: #possible that in that file there exist no line for that low_iso_abundance isotopologue
            insert_hitemp(iso_filepaths[i], isotop_names[i], line_source_id, hitemp_ref_link)           
            #delete the file because i am using appending method to open(file) and multiple files exists for one molecule
            os.remove(iso_filepaths[i])
            
        
    print("Finished id %s with %s isotopologues in %s seconds" % (M_id, len(isotop_names), time.time() - t))
    
##########################
    
#cant deal with ftp: manually downlaoding all files

start_time = time.time()

hitemp_ref_link = r'https://www.ucl.ac.uk/amopp/sites/amopp/files/480.pdf'
version_name_1 = 'HITEMP_2010'
version_name_2 = 'HITEMP_2019'

'''
#used to unzip all the downloaded files
for fp in os.listdir('/home/toma/Downloads'): 
    if fp.endswith('.zip'):
        with zipfile.ZipFile('/home/toma/Downloads/' + fp, 'r') as zip_ref:
            zip_ref.extractall('/home/toma/Desktop/linelists-database/hitemp')    
'''

hitemp_versions = ['HITEMP_2010', 'HITEMP_2019']
hitemp_ids = []
for version in hitemp_versions: 
    #insert the line_source into source_properties and get line_source_id
    insert_version_query = "INSERT INTO source_properties(line_source, max_temperature, max_nu, num_lines, bool_air, \
    bool_H2, bool_He, reference_link, line_source_id) VALUES('%s', null, null, null, 'YES', 'NO', 'NO', '%s', null);" % \
    (version, hitemp_ref_link)
    
    sql_order(insert_version_query)
    
    get_line_source_id_query = "SELECT line_source_id FROM source_properties WHERE line_source = '{}'".format(version)
    
    data = fetch(get_line_source_id_query)
    
    if len(data) != 1:
        raise Exception('should have exactly one line_source_id corresponding to one line_source')
        
    hitemp_ids.append(data[0][0])

#load molecule properties file
mol_ids, iso_ids, iso_names, iso_abundances, iso_masses, mol_names = \
np.loadtxt('/home/toma/Desktop/molecule_properties (copy).txt', dtype='str', skiprows=1, usecols=(1, 2, 3, 4, 5, 6), unpack=True)

counter = 0
for hitemp_fp in os.listdir('/home/toma/Desktop/linelists-database/hitemp'):
    counter += 1
    print(hitemp_fp)
    insert_all('/home/toma/Desktop/linelists-database/hitemp/' + hitemp_fp, hitemp_fp)
    print('Finished file', str(counter))

print("Finished in %s seconds" % (time.time() - start_time))