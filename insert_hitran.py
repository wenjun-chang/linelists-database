#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 13:13:51 2019

@author: toma
"""

#this code helps create tables and import hitran data into the database
#(did not do drop database or use database...add later)

import MySQLdb
from query_functions import sql_order

###############

def create_database():
    db = MySQLdb.connect(host='localhost', user='toma', passwd='Happy810@')
    cursor = db.cursor()
    cursor.execute('drop database if exists linelist')
    cursor.execute('create database linelist')
    db.commit()
    cursor.close()
    db.close()
        
#molecule_name format for example, CO2, is (13C)(16O)2
particles_table = '''create table if not exists particles (molecule_name varchar(12) not null, \
iso_name varchar(20) not null, iso_abundance double not null, iso_mass double not null, particle_id \
int unsigned not null auto_increment primary key); '''

#create table 2 for all the lines for each particle in table 1
#nu stands for transition wavenumber
#a stands for einstein coefficient
#gp stands for the degeneracy of the lower state
#is gp good with smallint??? should H2 He stuff default as null or 0.0
lines_table = '''create table if not exists transitions (nu double not null, A double not null, \
gamma_air double null, n_air double null, delta_air double null, \
elower double not null, gp smallint not null, gamma_H2 double null, \
n_H2 double null, delta_H2 double null, gamma_He double null, n_He double null, \
delta_He float null, data_type enum('HITRAN', 'EXOMOL') not null, \
version varchar(10) not null, particle_id int not null, line_id int unsigned \
not null auto_increment primary key); '''

#create table 3 for the partition coefficient across all temperatures for each
#particle in table 1
#should temoerature in K be float or int? also partition
partitions_table = '''create table if not exists partitions (temperature float not null, `partition` float not null, \
particle_id int not null, partition_id int unsigned not null auto_increment \
primary key); '''

#insert CO data into table particle
CO = '''insert ignore into particles values ('%s', '%s', '%s', '%s', null); ''' % ('CO', '(12C)(16O)', 0.986544, 27.994915)
#why the ignore statement is not working????
##########################

def insert_hitran(filename): 
    #connect to the database
    db = MySQLdb.connect(host='localhost', user='toma', passwd='Happy810@', db='linelist')
    #do put actual password when run
    
    #create a cursor object
    cursor = db.cursor()
    
    #insert the data of all lines for CO into table lines
    # with open('CO(copy).out') as infile: #
    try:
        #open the file
        infile = open(filename)
        for line in infile:
            data = line.split()
            for i in range(len(data)): 
                if data[i] == '#':
                    data[i] = None

            query = "insert into transitions values(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, 'HITRAN', '2012', 1, null)"
            cursor.execute(query, data)
        #maybe print the id each time to mnake sure it runs correctly?
        #should He and H2 stuff default as null or 0.0?!!!!!!!!!!
        #line table arrangement corresponding to tuple indexes: 
        #(nu, a, gamma_air, n_air, delta_air, elower, gp, gamma_H2, n_H2, delta_H2, gamma_He, n_He, delta_He, data_type, version, particle_id, line_id)
        #( 0, 1,      2,       3,       4,        5,    6,    7,       8,      9,       10,      11,    12,        13,       14,        15,        16  )
        #commit changes
        db.commit()
        infile.close()
    
    except Exception as e:
        #if errors occur
        db.rollback()
        print('insert hitran data failed', e)
      
    finally: 
        #close up cursor and connection
        cursor.close()
        db.close()
        
#################
        
def insert_exomol():
    pass

####################
        
def insert_partition(filename):
    #connect to the database
    db = MySQLdb.connect(host='localhost', user='toma', passwd='Happy810@', db='linelist')
    #do put actual password when run
    
    #create a cursor object
    cursor = db.cursor()
    
    #insert all the pratition functions into the table
    try:
        #open the file
        infile = open(filename)
        for line in infile:
           data = line.split()
           query = ('''insert into partitions (temperature, `partition`, particle_id, \
                                               partition_id) values('%s', '%s', 1, null)''' % (data[0], data[1]))
           cursor.execute(query) #maybe print the id each time to mnake sure it runs correctly?
        #commit changes
        db.commit()
        infile.close()
    except Exception as e: 
        #if errors occur
        db.rollback()
        print('insert partition failed', e)
    finally: 
        #close up cursor and connection
        cursor.close()
        db.close()

##################
        
def main():
    
    create_database()
    #create the tables
    #can i combine them? probably
    sql_order(particles_table)
    sql_order(lines_table)
    sql_order(partitions_table)
    
    #insert CO data in table 1
    sql_order(CO)
    
    #insert the data of all lines for CO into table lines
    insert_hitran('/home/toma/Desktop/co_test.out')
    
    #insert all the pratition functions into the table
    insert_partition('/home/toma/Desktop/12C-16O__Li2015_partition.pf')
    
    
    
    
    
if __name__ == '__main__':
    main()