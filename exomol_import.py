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
import itertools

DEFAULT_GAMMA = 0.0700
DEFAULT_N = 0.500

#################

start_time = time.time()

#connect to the database
db = MySQLdb.connect(host='localhost', user='toma', passwd='Happy810@', db='linelist') 

#create a cursor object
cursor = db.cursor()

#disable autocommit to improve performance

sql_order('SET autocommit = 0')
sql_order('SET unique_checks = 0')
sql_order('SET foreign_key_checks = 0')

#get parameters needed to insert exomol data into transitions

#states in id order starts in 1 
Es, gs, Js = np.loadtxt('/home/toma/Desktop/12C-16O__Li2015.states', usecols=(1, 2, 3), unpack=True)

#J starts with 0
gamma_H2s, n_H2s = np.loadtxt('/home/toma/Desktop/12C-16O__H2.broad', usecols=(1, 2), unpack=True)

#J starts with 0
gamma_Hes, n_Hes = np.loadtxt('/home/toma/Desktop/12C-16O__He.broad', usecols=(1, 2), unpack=True)


def insert_exomol(start_line, end_line, filename):
    upper_ids, lower_ids, As = np.loadtxt(itertools.islice(trans, start_line, end_line), usecols=(0, 1, 2), unpack=True)
    
    #file that the parameters are written into and import to mysql using LOAD DATA INFILE
    f = open(filename, 'w') 
    
    file_time = time.time()
    
    #need optimizeoptimize
    for i in range(len(upper_ids)):
        upper_id = int(upper_ids[i])
        lower_id = int(lower_ids[i])
        A = As[i]
        
        E_upper = Es[upper_id - 1]
        E_lower = Es[lower_id - 1]
        gp = gs[lower_id - 1]
        J_lower = int(Js[lower_id - 1])
        
        #store length into variable to run faster
        len_H2 = len(n_H2s)
        len_He = len(n_Hes)
    
        if J_lower < min(len_H2, len_He):
            gamma_H2 = gamma_H2s[J_lower]
            n_H2 = n_H2s[J_lower]
            gamma_He = gamma_Hes[J_lower]
            n_He = n_Hes[J_lower]
        
        elif len_H2 >= len_He and (J_lower >= len_He and J_lower < len_H2):
                gamma_H2 = gamma_H2s[J_lower]
                n_H2 = n_H2s[J_lower]
                gamma_He = DEFAULT_GAMMA
                n_He = DEFAULT_N
                
        elif len_H2 < len_He and (J_lower >= len_H2 and J_lower < len_He):
                gamma_H2 = DEFAULT_GAMMA
                n_H2 = DEFAULT_N
                gamma_He = gamma_Hes[J_lower]
                n_He = n_Hes[J_lower]
        else: 
            gamma_H2 = DEFAULT_GAMMA
            n_H2 = DEFAULT_N
            gamma_He = DEFAULT_GAMMA
            n_He = DEFAULT_N
        
        '''
        else:
            if len_H2 >= len_He:
                if J_lower >= len_He and J_lower < len_H2:
                    gamma_H2 = gamma_H2s[J_lower]
                    n_H2 = n_H2s[J_lower]
                    gamma_He = DEFAULT_GAMMA
                    n_He = DEFAULT_N
                else: 
                    gamma_H2 = DEFAULT_GAMMA
                    n_H2 = DEFAULT_N
                    gamma_He = DEFAULT_GAMMA
                    n_He = DEFAULT_N
            elif len_H2 < len_He:
                if J_lower >= len_H2 and J_lower < len_H2:
                    gamma_H2 = DEFAULT_GAMMA
                    n_H2 = DEFAULT_N
                    gamma_He = gamma_Hes[J_lower]
                    n_He = n_Hes[J_lower]
                else: 
                    gamma_H2 = DEFAULT_GAMMA
                    n_H2 = DEFAULT_N
                    gamma_He = DEFAULT_GAMMA
                    n_He = DEFAULT_N
        '''

        if gamma_H2 is None: 
            gamma_H2 = DEFAULT_GAMMA
        if gamma_He is None: 
            gamma_He = DEFAULT_GAMMA
        if n_H2 is None: 
            n_H2 = DEFAULT_N
        if n_He is None: 
            n_He = DEFAULT_N
        
        v_ij = E_upper - E_lower
        
        data = [v_ij, A, E_lower, gp, gamma_H2, n_H2, gamma_He, n_He]
            
        for item in data: 
            f.write("%s " % item)
        f.write("\n")
        
        global counter
        counter += 1
        #print(counter)
    
    f.close()
    
    print("Write infile in %s seconds" % (time.time() - file_time))

    load_time = time.time()
    
    print("Bulk inserting exomol data...")
    
    cursor.execute("LOAD DATA LOCAL INFILE '/home/toma/Desktop/exomol.txt' INTO TABLE transitions FIELDS TERMINATED BY ' ' LINES TERMINATED BY '\n' \
              (@col1, @col2, @col3, @col4, @col5, @col6, @col7, @col8) SET nu=@col1, A=@col2, elower=@col3, gp=@col4, \
              gamma_H2=@col5, n_H2=@col6, gamma_He=@col7, n_He=@col8, line_source='EXOMOL_Li2015', particle_id=1;")
    
    db.commit()
    
    print('Executed {} lines of exomol data'.format(counter))
    print("Bulk inserted exomol data in %s seconds" % (time.time() - load_time))    

#transition file
#max_rows =.....
#get the number of lines in trans file
length_trans = sum(1 for line in open('/home/toma/Desktop/12C-16O__Li2015.trans'))
print(length_trans)
counter = 0 
with open('/home/toma/Desktop/12C-16O__Li2015.trans') as trans:
    start_line = 0
    max_size = 1e5 
    while length_trans >= start_line + max_size:
        
        insert_exomol(start_line, int(start_line + max_size), '/home/toma/Desktop/exomol.txt')
        
        #start_line += 1e5
        #...islice removes all the lines used from the file
        length_trans -= max_size
        print(length_trans)
    

    #out of the while loop when difference between start_line and the max lines in trans file is less than 1e6
    insert_exomol(start_line, int(length_trans), '/home/toma/Desktop/exomol.txt')

trans.close()

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

print("Bulk inserting partition data...")
print("Executed {} lines of partition data".format(counter))

sql_bulk_order(query_insert_partitions, partition_data)    
db.commit()

################

#turn them back on

sql_order('SET unique_checks = 1')
sql_order('SET foreign_key_checks = 1')

cursor.close()
db.close()

print("Finished in %s seconds" % (time.time() - start_time))