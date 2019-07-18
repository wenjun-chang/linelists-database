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
CO = "INSERT INTO particles VALUES('%s', '%s', '%s', '%s', '%s', null);" % ('CO', '(12C)(16O)', 0.986544, 27.994915, 'HITRAN_2016')

##########################

def insert_hitran(filename, version_name, particle_id): 
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
            #(nu, a, gamma_air, n_air, delta_air, elower, g_upper, gamma_H2, n_H2, delta_H2, gamma_He, n_He, delta_He, line_source, particle_id, line_id)
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
                  delta_He=@col13, line_source= {}, particle_id={};".format(version_name, particle_id))
        
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
        
##################        
def main():
    
    start_time = time.time()
       
    #insert CO data in table 1
    #sql_order(CO)
    #PH3 = "INSERT INTO particles VALUES('%s', '%s', '%s', '%s', '%s', null);" % ('PH3', '(31P)(1H)3', 0.999533, 33.997238, 'EXOMOL_SAlTY')
    #sql_order(PH3)
    #insert the data of all lines for CO into table lines
    insert_hitran('/home/toma/Desktop/co_test.out', 'HITRAN_2016', 1)
        
    print("Finished in %s seconds" % (time.time() - start_time))
    
if __name__ == '__main__':
    main()