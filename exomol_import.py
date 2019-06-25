#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 16:02:30 2019

@author: toma
"""


#this code helps create tables needed to import exomol data and import exomol data

import MySQLdb
import numpy as np
from query_functions import sql_order

DEFAULT_GAMMA = 0.0700
DEFAULT_N = 0.500

#create table states needed to store the states file of exomol data for each molecule
#state_id is the exomol state_id, and id is the mysql primary key. 
#what should the size of state_id be??? smallint or mediumint????????
states_table_create_query = '''create table if not exists states (state_id mediumint not null, E double not null, \
g smallint not null, J smallint not null, particle_id int not null, id int unsigned \
not null auto_increment primary key); '''

#create table broad_params to store gamma, n, and delta for exomol data
broad_params_table_create_query =  '''create table if not exists broad_params \
(J smallint not null, gamma_H2 double null, n_H2 double null, \
gamma_He double null, n_He double null, particle_id int not null, \
broad_id int unsigned not null auto_increment primary key); '''

#create table 4 and 5
sql_order(states_table_create_query)
sql_order(broad_params_table_create_query)

#################

#connect to the database
db = MySQLdb.connect(host='localhost', user='toma', passwd='Happy810@', db='linelist') 
#do put actual password when run


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
        
#insert states file
insert_exomol_file(db, '/home/toma/Desktop/12C-16O__Li2015.states', \
                   "insert into states (state_id, E, g, J, particle_id, id) values({}, {}, {}, {}, 1 , null)",\
                   [0, 1, 2, 3])

#insert broad H2 file
insert_exomol_file(db, '/home/toma/Desktop/12C-16O__H2.broad', \
                   "insert into broad_params (J, gamma_H2, n_H2, gamma_He, n_He, \
                                          particle_id, broad_id) values({}, {}, {}, null, null, 1, null)", [3, 1, 2])

#insert broad He file
insert_exomol_file(db, '/home/toma/Desktop/12C-16O__He.broad', \
                   "update broad_params set gamma_He = {}, n_He = {} where J = {}", [1, 2, 3])


#####
#open trans file and fetch the data in table 4 and 5 according to state_id's in trans file
trans = open('/home/toma/Desktop/12C-16O__Li2015.trans')
counter = 0

#create a cursor object
cursor = db.cursor()

for line in trans:
    data = line.split()
        
    A = data[2]
    #get E_higher using higher_state_id
    upper = "select E from states where state_id = {}".format(data[0])
    cursor.execute(upper)
    E_higher, = cursor.fetchone()
        
    #get E_lower, gp, J using higher_state_id
    lower = "select E, g, J from states where state_id = {}".format(data[1])
    cursor.execute(lower)
    E_lower, gp, J = cursor.fetchone()
        
    v_ij = E_higher - E_lower
        #get other info needed from table broad_params using J
    get_param = "select gamma_H2, n_H2, gamma_He, n_He, particle_id from broad_params where J = {}".format(J)
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
    query = ("insert into transitions (nu, A, gamma_air, n_air, delta_air, elower, gp, gamma_H2, n_H2, \
                                         delta_H2, gamma_He, n_He, delta_He, data_type, version, particle_id, \
                                         line_id) values('%s', '%s', null, null, null, '%s', '%s', '%s', '%s', 0, '%s', '%s', 0, '%s', '%s', '%s', null)" \
                                        % (v_ij, A, E_lower, gp, param[0], param[1], param[2], param[3], 'EXOMOL', 'Li2015', param[4])) 
    #print(query)
    cursor.execute(query)
    counter += 1
    print("Processing line", counter, param)
#commit changes
db.commit()
trans.close()
     
cursor.close()
db.close()














