#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 16:02:30 2019

@author: toma
"""


#this code helps import exomol data

import MySQLdb
import numpy as np
from query_functions import sql_bulk_order, sql_order
import time

DEFAULT_GAMMA = 0.0700
DEFAULT_N = 0.500

#################

start_time = time.time()

#connect to the database
db = MySQLdb.connect(host='localhost', user='toma', passwd='Happy810@', db='linelist') 

#create a cursor object
cursor = db.cursor()


'''LOAD DATA LOCAL INFILE '/home/toma/Desktop/hitran.txt' INTO TABLE transitions FIELDS TERMINATED BY '\b' LINES TERMINATED BY '\r\n';'''
'''LOAD DATA LOCAL INFILE '/home/toma/Desktop/exomol.txt' INTO TABLE transitions FIELDS TERMINATED BY '\b' LINES TERMINATED BY '\r\n';'''

#disable autocommit to improve performance
sql_order('SET autocommit = 0')

#get parameters needed to insert exomol data into transitions

#states in id order starts in 1 
Es, gs, Js = np.loadtxt('/home/toma/Desktop/12C-16O__Li2015.states', usecols=(1, 2, 3), unpack=True)

#J starts with 0
gamma_H2s, n_H2s = np.loadtxt('/home/toma/Desktop/12C-16O__H2.broad', usecols=(1, 2), unpack=True)

#J starts with 0
gamma_Hes, n_Hes = np.loadtxt('/home/toma/Desktop/12C-16O__He.broad', usecols=(1, 2), unpack=True)

#transition file
upper_ids, lower_ids, As = np.loadtxt('/home/toma/Desktop/12C-16O__Li2015.trans', usecols=(0, 1, 2), unpack=True)

exomol_data = []
query_insert_exomol = "INSERT INTO transitions (nu, A, gamma_air, n_air, delta_air, elower, gp, \
gamma_H2, n_H2, delta_H2, gamma_He, n_He, delta_He, line_source, particle_id, line_id) \
VALUES(%s, %s, null, null, null, %s, %s, %s, %s, null, %s, %s, null, %s, %s, null)"

counter = 0

#
f = open('/home/toma/Desktop/exomol.txt', 'w') 

for i in range(len(upper_ids)):
    upper_id = int(upper_ids[i])
    lower_id = int(lower_ids[i])
    A = As[i]
    
    E_upper = Es[upper_id - 1]
    E_lower = Es[lower_id - 1]
    gp = gs[lower_id - 1]
    J_lower = int(Js[lower_id - 1])
    
    if J_lower < min(len(n_H2s), len(n_Hes)):
        gamma_H2 = gamma_H2s[J_lower]
        n_H2 = n_H2s[J_lower]
        gamma_He = gamma_Hes[J_lower]
        n_He = n_Hes[J_lower]
        
    else:
        if len(n_H2s) >= len(n_Hes):
            if J_lower >= len(n_Hes) and J_lower < len(n_H2s):
                gamma_H2 = gamma_H2s[J_lower]
                n_H2 = n_H2s[J_lower]
                gamma_He = DEFAULT_GAMMA
                n_He = DEFAULT_N
            else: 
                gamma_H2 = DEFAULT_GAMMA
                n_H2 = DEFAULT_N
                gamma_He = DEFAULT_GAMMA
                n_He = DEFAULT_N
        elif len(n_H2s) < len(n_Hes):
            if J_lower >= len(n_H2s) and J_lower < len(n_Hes):
                gamma_H2 = DEFAULT_GAMMA
                n_H2 = DEFAULT_N
                gamma_He = gamma_Hes[J_lower]
                n_He = n_Hes[J_lower]
            else: 
                gamma_H2 = DEFAULT_GAMMA
                n_H2 = DEFAULT_N
                gamma_He = DEFAULT_GAMMA
                n_He = DEFAULT_N
            
    if gamma_H2 is None: 
        gamma_H2 = DEFAULT_GAMMA
    if gamma_He is None: 
        gamma_He = DEFAULT_GAMMA
    if n_H2 is None: 
        n_H2 = DEFAULT_N
    if n_He is None: 
        n_He = DEFAULT_N
    
    v_ij = E_upper - E_lower
    
    #
    data = [v_ij, A, 'null', 'null', 'null', E_lower, gp, gamma_H2, n_H2, 'null', gamma_He, n_He, 'null', 'EXOMOL_Li2015', 1, 'null']
    for item in data: 
        f.write("%s " % item)
    f.write("\n")
    
    
    '''
    exomol_data.append((v_ij, A, E_lower, gp, gamma_H2, n_H2, gamma_He, n_He, 'EXOMOL_Li2015', 1))
    
    
    
    counter += 1
    #print("Processing line {} for exomol data".format(counter))

print("Bulk inserting exomol data...")

bulk_time = time.time()

#test speed for mega data

for i in range(13):
    exomol_data += exomol_data
    print(len(exomol_data))


sql_bulk_order(query_insert_exomol, exomol_data)
db.commit()

print("Bulk inserted exomol data in %s seconds" % (time.time() - bulk_time))

##################

#insert partition file
Ts, partition_functions = np.loadtxt('/home/toma/Desktop/12C-16O__Li2015_partition.pf', usecols=(0, 1), unpack=True)

partition_data = [] 
query_insert_partitions = "INSERT INTO partitions (temperature, `partition`, particle_id, partition_id) VALUES(%s, %s, 1, null)"

counter = 0
for j in range(len(partition_functions)):
    T = Ts[j]
    partition = partition_functions[j]
    
    partition_data.append((T, partition))
    
    counter += 1
    #print("Processing line {} for partition data".format(counter))

print("Bulk inserting partition data...")
sql_bulk_order(query_insert_partitions, partition_data)    
db.commit()

################

cursor.close()
db.close()

print("Finished in %s seconds" % (time.time() - start_time))

'''