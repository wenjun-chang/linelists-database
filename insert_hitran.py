#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 13:13:51 2019

@author: toma
"""

#this code helps insert hitran data into the database

import MySQLdb
from query_functions import sql_order, fetch

##########################

def insert_hitran(filename, version_name, particle_id, reference_link): 
    #connect to the database
    db = MySQLdb.connect(host='localhost', user='toma', passwd='Happy810@', db='linelist')
    
    #create a cursor object
    cursor = db.cursor()
    
    #disable autocommit to improve performance
    sql_order('SET autocommit = 0')
    sql_order('SET unique_checks = 0')
    sql_order('SET foreign_key_checks = 0')
    sql_order('SET sql_log_bin = 0')
    
    #insert the data of all lines for CO into table lines
    # with open('CO(copy).out') as infile: #
    try:
        #insert the line_source into source_properties and get line_source_id
        insert_version_query = "INSERT IGNORE INTO source_properties(line_source, max_temperature, max_nu, num_lines, bool_air, \
        bool_H2, bool_He, reference_link, particle_id, line_source_id) VALUES('%s', null, null, null, 'YES', 'YES', 'YES', '%s', \
        '%s', null);" % (version_name, reference_link, particle_id)
            
        sql_order(insert_version_query)
        
        get_line_source_id_query = "SELECT line_source_id FROM source_properties WHERE line_source = '{}' AND \
        particle_id = {}".format(version_name, particle_id)
        
        data = fetch(get_line_source_id_query)
        
        if len(data) != 1:
            raise Exception('should have exactly one line_source_id corresponding to one line_source')
            
        line_source_id = data[0][0]

        #file that the parameters are written into and import to mysql using LOAD DATA INFILE
        f = open('/home/toma/Desktop/hitran.txt', 'w') 
        
        #open the file
        infile = open(filename)
        
        counter = 0
        
        #create a list of all the queries to bulk insert it
        #bulk_data = []
        #query = "INSERT INTO transitions VALUES(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, 'HITRAN_2016', 1, null)"
        
        for line in infile:
            data = line.strip().split(',')
            for i in range(len(data)): 
                if data[i] == '#':
                    data[i] = '\\N'
            #line table arrangement corresponding to tuple indexes: 
            #(nu, a, gamma_air, n_air, delta_air, elower, g_upper, gamma_H2, n_H2, delta_H2, gamma_He, n_He, delta_He, line_source_id, particle_id, line_id)
            #( 0, 1,      2,       3,       4,        5,    6,    7,       8,      9,       10,      11,    12,        13,         14,         15  )
            
            #make sure at least one gamma and one n value is not null
            if data[2] == '\\N' and data[7] == '\\N' and data[10] == '\\N':
                raise Exception('should have at least one gamma value')
            if data[3] == '\\N' and data[8] == '\\N' and data[11] == '\\N':
                raise Exception('should have at least one n value')
            
            
            #write into infile the parameters for each line
            for item in data: 
                f.write("%s " % item)
            f.write("\n")
            
            counter += 1
        
        f.close()
        print("Bulk inserting hitran data...")
        
        cursor.execute("LOAD DATA LOCAL INFILE '/home/toma/Desktop/hitran.txt' INTO TABLE transitions FIELDS TERMINATED BY ' ' LINES TERMINATED BY '\n' \
                  (@col1, @col2, @col3, @col4, @col5, @col6, @col7, @col8, @col9, @col10, @col11, @col12, @col13) SET nu=@col1, A=@col2, gamma_air=@col3, \
                  n_air=@col4, delta_air=@col5, elower=@col6, g_upper=@col7, gamma_H2=@col8, n_H2=@col9, delta_H2=@col10, gamma_He=@col11, n_He=@col12, \
                  delta_He=@col13, line_source_id={}, particle_id={};".format(line_source_id, particle_id))
        
        #commit changes and close file
        db.commit()
        infile.close()
        
        #turn it back on
        sql_order('SET unique_checks = 1')
        sql_order('SET foreign_key_checks = 1')
        sql_order('SET sql_log_bin = 1')
        
        print('Executed {} lines of hitran data'.format(counter))
        
    except Exception as e:
        #if errors occur
        db.rollback()
        print('insert hitran data failed', e)
      
    finally: 
        #close up cursor and connection
        cursor.close()
        db.close()
        