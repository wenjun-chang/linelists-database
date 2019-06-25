#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 16:02:30 2019

@author: toma
"""


#this code helps create tables needed to import exomol data and import exomol data

import MySQLdb

######################

def sql_order(query):
    #connect to the database
    db = MySQLdb.connect(host='localhost', user='toma', passwd='Happy810@', db='linelist')
    #do put actual password when run
    
    #create a cursor object
    cursor = db.cursor()
    #create table 1 for particles
    #create the tables
    try: 
        #execute order in mysql
        cursor.execute(query)
        #commit changes
        db.commit()
    except: 
        #if errors occur
        db.rollback()
        
    finally: 
        #close up cursor and connections
        cursor.close()
        db.close()

################
        
#create table states needed to store the states file of exomol data for each molecule
#state_id is the exomol state_id, and id is the mysql primary key. 
#what should the size of state_id be??? smallint or mediumint????????
table4 = '''create table if not exists `states` (`state_id` mediumint not null, `E` double not null, \
`g` smallint not null, `J` smallint not null, particle_id int not null, `id` int unsigned \
not null auto_increment primary key); '''

#create table broad_param to store gamma, n, and delta for exomol data
table5 =  '''create table if not exists broad_param (J smallint not null, gamma_H2 double not null default \
0.0, n_H2 double not null default 0.0, delta_H2 double not null default 0.0, \
gamma_He double not null default 0.0, n_He double not null default 0.0, \
delta_He double not null default 0.0, particle_id int not null, broad_id int unsigned \
not null auto_increment primary key); '''

#create table 4 and 5
sql_order(table4)
sql_order(table5)

#################

#connect to the database
db = MySQLdb.connect(host='localhost', user='toma', passwd='Happy810@', db='linelist') 
#do put actual password when run

#create a cursor object
cursor = db.cursor()
'''
try:
'''
#insert states file
infile = open('/home/toma/Desktop/12C-16O__Li2015.states')
for line in infile:
    data = line.split(' ')
    query = ('''insert into `states` (state_id mediumint not null, `E` double not null, \
                                        `g` smallint not null, `J` smallint not null, particle_id int not null, `id` int unsigned \
                                        not null auto_increment primary key) values('%s', '%s', '%s', '%s', '%s', '%s')''' \
                                        % (data[0], data[1], data[2], data[3], 1, 'null')) 

    cursor.execute(query)
    #maybe print the id each time to mnake sure it runs correctly?
    #commit changes
    db.commit()
infile.close()
    
#insert broad H2 file
infile2 = open('/home/toma/Desktop/12C-16O__H2.broad')
for line in infile2:
    data = line.split(' ')
    query = ('''insert into broad_param (J smallint not null, gamma_H2 double null default \
                                             0.0, n_H2 double null default 0.0, delta_H2 double null default 0.0, \
                                             gamma_He double null default 0.0, n_He double null default 0.0, \
                                             delta_He double null default 0.0, particle_id int not null, broad_id int unsigned \
                                             not null auto_increment primary key) values('%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s')''' \
                                            % (data[3], data[1], data[2], 'null', 'null', 'null', 'null' , 1, 'null')) 

    cursor.execute(query)
    #maybe print the id each time to mnake sure it runs correctly?
    #commit changes
db.commit()
infile2.close()
    
#insert broad He file
infile3 = open('/home/toma/Desktop/12C-16O__He.broad')
for line in infile3:
    data = line.split(' ')
    query = ('''update broad_param set gamma_He = {}, n_He double = {} where J = {}'''.format(data[1], data[2], data[3]))

    cursor.execute(query)
    #maybe print the id each time to mnake sure it runs correctly?
    #commit changes
db.commit()
infile3.close()
    
#insert delta H2 file
#not sure if this file is for delta and if delta = beta in the file?!!!
infile4 = open('/home/toma/Desktop/fit-co.tab')
for line in infile4:
    data = line.split(' ')
    query = ('''update broad_param set delta_H2 = {} where J = {}'''.format(data[4], data[1]))

    cursor.execute(query)
    #maybe print the id each time to mnake sure it runs correctly?
    #commit changes
db.commit()
infile4.close()
    
    
    #####
    #open trans file and fetch the data in table 4 and 5 according to state_id's in trans file
trans = open('/home/toma/Desktop/12C-16O__Li2015.trans')
for line in trans:
    data = line.split(' ')
        
    A = data[2]
    #get E_higher using higher_state_id
    upper = '''select E from states where state_id = {}'''.format(data[0])
    cursor.execute(query)
    temp = cursor.fetchone()
    E_higher = temp[0]
        
    #get E_lower, gp, J using higher_state_id
    lower = '''select E, g, J from states where state_id = {}'''.format(data[1])
    cursor.execute(query)
    temp2 = cursor.fetchone()
    E_lower = temp2[0]
    gp = temp2[1]
    J = temp2[2]
        
    v_ij = E_lower - E_higher
        #get other info needed from table broad_param using J
    get_param = '''select * from broad_param where J = {}'''.format(J)
    cursor.execute(query)
    param = cursor.fetchone()
        #(J, gamma_H2, n_H2, delta_H2, gamma_He, n_He, delta_He, particle_id, param_id)
        #(0,   1,       2,      3,       4,       5,     6,            7,        8    )
        
        #put all the exomol data into table lines
    query = ('''insert into lines (nu double not null, a double not null, \
                                       gammar_air double null, n_air double null, delta_air double null, \
                                       elower double not null, gp smallint not null, gamma_H2 double not null default \
                                       0.0, n_H2 double not null default 0.0, delta_H2 double not null default 0.0, \
                                       gamma_He double not null default 0.0, n_He double not null default 0.0, \
                                       delta_He double not null default 0.0, data_type enum('HITRAN', 'EXOMOL') not null, \
                                       version varchar(10) not null, particle_id int not null, line_id int unsigned \
                                       not null auto_increment primary key) values('%s', '%s', '%s', '%s', '%s',\
                                       '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s')''' \
                                        % (v_ij, A, 'null', 'null', 'null', E_lower, gp, param[1], param[2], param[3], param[4], \
                                           param[5], param[6], 'EXOMOL', 'Li2015', param[7], 'null')) 
        
        
#commit changes
db.commit()
trans.close()
'''
except: 
    #if errors occur
    db.rollback()

finally: 
        #close up cursor and connection'''
cursor.close()
db.close()















