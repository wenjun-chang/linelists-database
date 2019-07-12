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

################
        
#molecule_name format for example, CO2, is (13C)(16O)2
particles_table_create_query = "CREATE TABLE IF NOT EXISTS particles (molecule_name VARCHAR(5) NOT NULL, \
iso_name VARCHAR(25) NOT NULL, iso_abundance DOUBLE NOT NULL, iso_mass DOUBLE NOT NULL, default_line_source \
VARCHAR(25) NOT NULL, particle_id INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY); "

#create table for all the lines for each particle in table particles
#nu stands for transition wavenumber
#a stands for einstein coefficient
#g_upper stands for the degeneracy of the uppper state
transitions_table_create_query = "CREATE TABLE IF NOT EXISTS transitions (nu DOUBLE NOT NULL, A FLOAT NOT NULL, \
gamma_air FLOAT, n_air FLOAT, delta_air FLOAT, elower DOUBLE NOT NULL, g_upper SMALLINT NOT NULL, gamma_H2 FLOAT, \
n_H2 FLOAT, delta_H2 FLOAT, gamma_He FLOAT, n_He FLOAT, delta_He FLOAT, line_source VARCHAR(25) NOT NULL, \
particle_id INT UNSIGNED NOT NULL, line_id INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY, FOREIGN KEY transitions(particle_id) \
REFERENCES particles(particle_id) ON UPDATE CASCADE ON DELETE CASCADE) ROW_FORMAT=COMPRESSED;"

#create table for the partition coefficient across all temperatures for each particle in table 1
partitions_table_create_query = "CREATE TABLE IF NOT EXISTS partitions (temperature FLOAT NOT NULL, `partition` FLOAT NOT NULL, \
particle_id INT UNSIGNED NOT NULL, partition_id INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY, FOREIGN KEY partitions(particle_id) \
REFERENCES particles(particle_id) ON UPDATE CASCADE ON DELETE CASCADE);" 

#create table states needed to store the states file of exomol data for each molecule
#state_id is the exomol state_id, and id is the mysql primary key. 
#what should the size of state_id be??? smallint or mediumint????????
states_table_create_query = "CREATE TABLE IF NOT EXISTS states (state_id MEDIUMINT NOT NULL, E DOUBLE NOT NULL, \
g SMALLINT NOT NULL, J SMALLINT NOT NULL, particle_id INT UNSIGNED NOT NULL, id INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY);"

#create table broad_params to store gamma, n, and delta for exomol data
broad_params_table_create_query =  "CREATE TABLE IF NOT EXISTS broad_params \
(J SMALLINT NOT NULL, gamma_H2 DOUBLE, n_H2 DOUBLE, gamma_He DOUBLE, n_He DOUBLE, \
particle_id INT UNSIGNED NOT NULL, broad_id INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY);"

##################

#create indexes on nu, A, elower, and line_source in table transitions
create_index_line_source = "CREATE INDEX line_source_index ON transitions(line_source) USING BTREE;"
create_index_nu = "CREATE INDEX nu_index ON transitions(nu) USING BTREE;"
create_index_A = "CREATE INDEX A_index ON transitions(A) USING BTREE;"
create_index_elower = "CREATE INDEX elower_index ON transitions(elower) USING BTREE;"

##################
        
def main():
    
    start_time = time.time()
    
    #create the database first and drop it if it exists already
    create_database()
    
    #create the tables
    sql_order(particles_table_create_query)
    sql_order(transitions_table_create_query)
    sql_order(partitions_table_create_query)
    sql_order(states_table_create_query)
    sql_order(broad_params_table_create_query)
    
    
    '''
    #create the indexes in table transitions
    #index significance: line_source > nu > A > elower
    sql_order(create_index_line_source)
    sql_order(create_index_nu)
    sql_order(create_index_A)
    sql_order(create_index_elower)
    '''
   
    print("Finished in %s seconds" % (time.time() - start_time))
        
if __name__ == '__main__':
    main()