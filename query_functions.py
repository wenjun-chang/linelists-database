#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 11:56:36 2019

@author: toma
"""

import MySQLdb

##############

def sql_order(query):
    #connect to the database
    db = MySQLdb.connect(host='localhost', user='toma', passwd='Happy810@', db='linelist') 
    #do put actual password when run
    
    #create a cursor object
    cursor = db.cursor()
    
    try: 
        #execute order in mysql
        cursor.execute(query)
        #commit changes
        db.commit()
    except Exception as e: 
        #if errors occur
        db.rollback()
        print(e)
        
    finally: 
        #close up cursor and connections
        cursor.close()
        db.close()
    
################