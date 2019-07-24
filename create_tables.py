#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 10:49:14 2019

@author: toma
"""
#this file creates database and tables in mysql

import MySQLdb
from query_functions import sql_order
import time

###############

def create_database():
    db = MySQLdb.connect(host='localhost', user='toma', passwd='Happy810@')
    cursor = db.cursor()
    cursor.execute('DROP DATABASE IF EXISTS linelist')
    cursor.execute('CREATE DATABASE linelist')
    db.commit()
    cursor.close()
    db.close()

##################

#molecule_name format for example, CO2, is (13C)(16O)2
particles_table_create_query = "CREATE TABLE IF NOT EXISTS particles (\
molecule_name VARCHAR(5) NOT NULL, \
iso_name VARCHAR(25) NOT NULL, \
iso_abundance DOUBLE NOT NULL, \
iso_mass DOUBLE NOT NULL, \
default_line_source_id \
SMALLINT NOT NULL, \
particle_id SMALLINT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY\
);"

#create table for all the lines for each particle in table particles
#nu stands for transition wavenumber
#a stands for einstein coefficient
#g_upper stands for the degeneracy of the uppper state
transitions_table_create_query = \
"CREATE TABLE IF NOT EXISTS transitions (\
nu DOUBLE NOT NULL, A FLOAT NOT NULL, \
gamma_air FLOAT, \
n_air FLOAT, \
delta_air FLOAT, \
elower DOUBLE NOT NULL, \
g_upper SMALLINT NOT NULL, \
gamma_H2 FLOAT, \
n_H2 FLOAT, \
delta_H2 FLOAT, \
gamma_He FLOAT, \
n_He FLOAT, \
delta_He FLOAT, \
line_source_id SMALLINT UNSIGNED NOT NULL, \
particle_id SMALLINT UNSIGNED NOT NULL, \
line_id BIGINT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY, \
FOREIGN KEY (particle_id) REFERENCES particles(particle_id) ON UPDATE CASCADE ON DELETE CASCADE, \
FOREIGN KEY (line_source_id) REFERENCES source_properties(line_source_id) ON UPDATE CASCADE ON DELETE CASCADE) \
ROW_FORMAT=COMPRESSED;"

#create table for the partition coefficient across all temperatures for each particle in table 1
partitions_table_create_query = \
"CREATE TABLE IF NOT EXISTS partitions (\
temperature FLOAT NOT NULL, \
`partition` FLOAT NOT NULL, \
line_source_id SMALLINT UNSIGNED NOT NULL, \
particle_id SMALLINT UNSIGNED NOT NULL, \
partition_id INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY, \
FOREIGN KEY (particle_id) REFERENCES particles(particle_id) ON UPDATE CASCADE ON DELETE CASCADE, \
FOREIGN KEY (line_source_id) REFERENCES source_properties(line_source_id) ON UPDATE CASCADE ON DELETE CASCADE\
);" 

#create table source_properties to store the limits and availbility of the parameters for each source from HITRAN or EXOMOL
source_properties_table_create_query = "CREATE TABLE IF NOT EXISTS source_properties (\
line_source VARCHAR(25) NOT NULL, \
max_temperature SMALLINT, \
max_nu DOUBLE, \
num_lines BIGINT, \
bool_air ENUM('YES', 'NO'), \
bool_H2 ENUM('YES', 'NO'), \
bool_He ENUM('YES', 'NO'), \
reference_link VARCHAR(250) NOT NULL, \
particle_id  SMALLINT UNSIGNED NOT NULL, \
line_source_id SMALLINT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY, \
FOREIGN KEY (particle_id) REFERENCES particles(particle_id) ON UPDATE CASCADE ON DELETE CASCADE);"

##################

#create indexes on nu, A, elower, and line_source in table transitions
create_index_line_source_id = "CREATE INDEX line_source_index ON transitions(line_source_id) USING BTREE;"
create_index_nu = "CREATE INDEX nu_index ON transitions(nu) USING BTREE;"
create_index_A = "CREATE INDEX A_index ON transitions(A) USING BTREE;"
create_index_elower = "CREATE INDEX elower_index ON transitions(elower) USING BTREE;"

##################
        
def main():
    
    start_time = time.time()
    '''
    #create the database first and drop it if it exists already
    create_database()
    
    #create the tables
    sql_order(particles_table_create_query)
    sql_order(source_properties_table_create_query)
    sql_order(transitions_table_create_query)
    sql_order(partitions_table_create_query)
    
    '''
    #Finished line source id index in 239.55348300933838 seconds
    #Finished nu index in 407.1092050075531 seconds
    #Finished A index in 332.6726574897766 seconds
    #Finished elower index in 459.4686324596405 seconds
    #Finished in 1438.8043451309204 seconds
    #for total 137802916 lines
    
    #create the indexes in table transitions
    #index significance: line_source_id > nu > A > elower
    print('starting...')
    t1 = time.time()
    sql_order(create_index_line_source_id)
    print("Finished line source id index in %s seconds" % (time.time() - t1))
    t2 = time.time()
    sql_order(create_index_nu)
    print("Finished nu index in %s seconds" % (time.time() - t2))
    t3 = time.time()
    sql_order(create_index_A)
    print("Finished A index in %s seconds" % (time.time() - t3))
    t4 = time.time()
    sql_order(create_index_elower)
    print("Finished elower index in %s seconds" % (time.time() - t4))
    
   
    print("Finished in %s seconds" % (time.time() - start_time))
        
if __name__ == '__main__':
    main()