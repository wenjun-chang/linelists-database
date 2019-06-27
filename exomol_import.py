#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 16:02:30 2019

@author: toma
"""


#this code helps import exomol data

import MySQLdb
import numpy as np
from query_functions import sql_order
import time

DEFAULT_GAMMA = 0.0700
DEFAULT_N = 0.500

#################

start_time = time.time()

#connect to the database
db = MySQLdb.connect(host='localhost', user='toma', passwd='Happy810@', db='linelist') 

def insert_exomol_file(db, filename, query_format, use_cols):
    try: 
        cursor = db.cursor()
        infile = open(filename)
        for line in infile:
            data = np.array(line.split())
            query = query_format.format(*data[np.array(use_cols)])
            print(query)
            cursor.execute(query)
       
        db.commit()
        infile.close()
    
    except Exception as e:
        db.rollback()
        print(e)
        
#disable autocommit to improve performance
sql_order('SET autocommit = 0')

#insert partitions file
insert_exomol_file(db, '/home/toma/Desktop/12C-16O__Li2015_partition.pf', \
                   "INSERT INTO partitions (temperature, `partition`, particle_id, \
                   partition_id) VALUES({}, {}, 1, null)", [0, 1])
       
#insert states file
insert_exomol_file(db, '/home/toma/Desktop/12C-16O__Li2015.states', \
                   "INSERT INTO states (state_id, E, g, J, particle_id, id) VALUES({}, {}, {}, {}, 1 , null)",\
                   [0, 1, 2, 3])

#insert broad H2 file
insert_exomol_file(db, '/home/toma/Desktop/12C-16O__H2.broad', \
                   "INSERT INTO broad_params (J, gamma_H2, n_H2, gamma_He, n_He, \
                                          particle_id, broad_id) VALUES({}, {}, {}, null, null, 1, null)", [3, 1, 2])

#insert broad He file
insert_exomol_file(db, '/home/toma/Desktop/12C-16O__He.broad', \
                   "UPDATE broad_params SET gamma_He = {}, n_He = {} WHERE J = {}", [1, 2, 3])


#####
#open trans file and fetch the data in table 4 and 5 according to state_id's in trans file
trans = open('/home/toma/Desktop/12C-16O__Li2015.trans')
counter = 0

#create a cursor object
cursor = db.cursor()

bulk_data = []

for line in trans:
    data = line.split()
        
    A = data[2]
    #get E_higher using higher_state_id
    upper = "SELECT E FROM states WHERE state_id = {}".format(data[0])
    cursor.execute(upper)
    E_higher, = cursor.fetchone()
        
    #get E_lower, gp, J using higher_state_id
    lower = "SELECT E, g, J FROM states WHERE state_id = {}".format(data[1])
    cursor.execute(lower)
    E_lower, gp, J = cursor.fetchone()
        
    v_ij = E_higher - E_lower
        #get other info needed from table broad_params using J
    get_param = "SELECT gamma_H2, n_H2, gamma_He, n_He, particle_id FROM broad_params WHERE J = {}".format(J)
    cursor.execute(get_param)
    
    
    ########need modify later
    param_tuple = cursor.fetchone()
    if param_tuple is None: 
        param = (DEFAULT_GAMMA, DEFAULT_N, DEFAULT_GAMMA, DEFAULT_N, 1)
    else: 
        param = list(param_tuple)
         
        #####
        #need to get this default value from web for each isotope later
        if param[1] is None: 
            param[1] = DEFAULT_GAMMA
        if param[3] is None: 
            param[3] = DEFAULT_GAMMA
        if param[2] is None: 
            param[2] = DEFAULT_N
        if param[4] is None: 
            param[4] = DEFAULT_N
    
        #(gamma_H2, n_H2, gamma_He, n_He, particle_id)
        #(   0,      1,      2,      3,       4      )
    
     
    #put all the exomol data into table lines
    query = ("INSERT INTO transitions (nu, A, gamma_air, n_air, delta_air, elower, gp, gamma_H2, n_H2, \
                                         delta_H2, gamma_He, n_He, delta_He, line_source, particle_id, \
                                         line_id) VALUES('%s', '%s', null, null, null, '%s', '%s', '%s', '%s', 0, '%s', '%s', 0, '%s', '%s', null)" \
                                        % (v_ij, A, E_lower, gp, param[0], param[1], param[2], param[3], 'EXOMOL Li2015', param[4])) 
    
    cursor.execute(query)
   
    
    counter += 1
    print("Processing line", counter, param)


#commit changes
db.commit()
trans.close()
     
cursor.close()
db.close()

print("Finished in %s seconds" % (time.time() - start_time))