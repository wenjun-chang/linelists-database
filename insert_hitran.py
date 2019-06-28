#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 13:13:51 2019

@author: toma
"""

#this code imports hitran data into the database

import MySQLdb
from query_functions import sql_order
import time

###############

#insert CO data into table particle
CO = "INSERT IGNORE INTO particles VALUES('%s', '%s', '%s', '%s', '%s', null);" % ('CO', '(12C)(16O)', 0.986544, 27.994915, 'HITRAN 2016')

##########################

def insert_hitran(filename): 
    #connect to the database
    db = MySQLdb.connect(host='localhost', user='toma', passwd='Happy810@', db='linelist')
    
    #create a cursor object
    cursor = db.cursor()
    
    #insert the data of all lines for CO into table lines
    # with open('CO(copy).out') as infile: #
    try:
        #open the file
        infile = open(filename)
        
        #create a list of all the queries to bulk insert it
        bulk_data = []
        counter = 0
        
        query = "INSERT INTO transitions VALUES(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, 'HITRAN 2016', 1, null)"
        
        for line in infile:
            data = line.split()
            for i in range(len(data)): 
                if data[i] == '#':
                    data[i] = None
            '''
            if data[4] is None: 
                data[4] = 0.0
            if data[9] is None: 
                data[9] = 0.0
            if data[12] is None: 
                data[12] = 0.0
     
            query = "INSERT INTO transitions VALUES(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, 'HITRAN 2016', 1, null)"
            cursor.execute(query, data)
            '''
            
            #line table arrangement corresponding to tuple indexes: 
            #(nu, a, gamma_air, n_air, delta_air, elower, gp, gamma_H2, n_H2, delta_H2, gamma_He, n_He, delta_He, line_source, particle_id, line_id)
            #( 0, 1,      2,       3,       4,        5,    6,    7,       8,      9,       10,      11,    12,        13,         14,         15  )
            
            if data[2] is None and data[7] is None and data[10] is None:
                raise Exception('should have at least one gamma value')
            if data[3] is None and data[8] is None and data[11] is None:
                raise Exception('should have at least one n value')
                
            
            bulk_data.append(tuple(data))
            
            counter += 1
            print("Processing line {} for hitran data".format(counter))
            
        #print(np.array(bulk_data).shape())
        print("Bulk inserting hitran data...")
        cursor.executemany(query, bulk_data)
        
        
        #commit changes and close file
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

##################
        
def main():
    
    start_time = time.time()
    
    #disable autocommit to improve performance
    sql_order('SET autocommit = 0')
    
    #insert CO data in table 1
    sql_order(CO)
    
    #insert the data of all lines for CO into table lines
    insert_hitran('/home/toma/Desktop/co_test.out')
    
    print("Finished in %s seconds" % (time.time() - start_time))
    
if __name__ == '__main__':
    main()