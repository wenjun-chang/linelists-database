#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 16:02:30 2019

@author: toma
"""

#this code helps import exomol data

import MySQLdb
import numpy as np
from query_functions import sql_bulk_order, sql_order, fetch
import time
import itertools

'''
#for CO
DEFAULT_GAMMA = 0.0700 #750 PH3
DEFAULT_N = 0.500 #530
particle_id = 30
'''
#for PH3
DEFAULT_GAMMA = 0.0750
DEFAULT_N = 0.530
particle_id = 93 # for PH3

####################

#insert partition file
def insert_partitions(partitions_filepath, line_source_id, particle_id): 
    Ts, partition_functions = np.loadtxt(partitions_filepath, usecols=(0, 1), unpack=True)
    partition_data = [] 
    query_insert_partitions = "INSERT INTO partitions (temperature, `partition`, line_source_id, particle_id, \
    partition_id) VALUES(%s, %s, {}, {}, null)".format(line_source_id, particle_id)
    
    counter = 0
    for j in range(len(partition_functions)):
        T = Ts[j]
        partition = partition_functions[j]
        
        partition_data.append((T, partition))
        
        counter += 1
    
    print("Bulk inserting partition data...")
    print("Executed {} lines of partition data".format(counter))
    
    sql_bulk_order(query_insert_partitions, partition_data)    
    db.commit()

#####################

#helper function for temporarily storing the data of H2/He broadening parameters
def temp_broad_param_dict(infile):
    broad_param_dict = {}
    with open(infile) as f: 
        #store [gamma_H2/Hes, n_H2/Hes] as a list in the dictionary with key as J + '_' + K
        for line in f: 
            data = line.strip().split()
            #if the prefix code is a1 or c1
            if data[0].endswith('1'):
                gamma = data[1]
                n = data[2]
                J = data[3]
                K = data[4]
                broad_param = {J + '_' + K : [gamma, n]}
                broad_param_dict.update(broad_param)
            elif data[0].endswith('0'):
                gamma = data[1]
                n = data[2]
                J = data[3]
                broad_param = {J : [gamma, n]}
                broad_param_dict.update(broad_param)
    return broad_param_dict

####################

def insert_exomol(start_line, end_line, outfile_name, infile, line_source_id, particle_id, no_broadening_param=False):
    
    upper_ids, lower_ids, As = np.loadtxt(itertools.islice(infile, start_line, end_line), usecols=(0, 1, 2), unpack=True)
    print(len(upper_ids), 'lines : Loaded the parameters from the transition file')
    
    #file that the parameters are written into and import to mysql using LOAD DATA INFILE
    f = open(outfile_name, 'w') 
    
    file_time = time.time()
    
    #need optimizeoptimize
    for i in range(len(upper_ids)):
        upper_id = int(upper_ids[i])
        lower_id = int(lower_ids[i])
        A = As[i]
        
        E_upper = Es[upper_id - 1]
        E_lower = Es[lower_id - 1]
        g_upper = int(gs[upper_id - 1])
        
        if no_broadening_param is True: 
            gamma_H2 = DEFAULT_GAMMA
            n_H2 = DEFAULT_N
            gamma_He = DEFAULT_GAMMA
            n_He = DEFAULT_N
        
        else: #if H2/He broad_param files exists i.e. no_broadening_param is False
            J_lower = int(Js[lower_id - 1])
            #K for c1 or a1 param style
            K_lower = int(Ks[lower_id - 1])
            
            #get H2 params
            if H2_dict.get(str(J_lower) + '_' + str(K_lower)) is not None: #look for specific params by J and K
                H2_params = H2_dict.get(str(J_lower) + '_' + str(K_lower))
                if H2_params[0] is not None: 
                    gamma_H2 = H2_params[0]
                else: 
                    gamma_H2 = DEFAULT_GAMMA
                if H2_params[1] is not None: 
                    n_H2 = H2_params[1]
                else: 
                    n_H2 = DEFAULT_N
            elif H2_dict.get(str(J_lower)) is not None: #use the more generalized params by J
                H2_params = H2_dict.get(str(J_lower))
                if H2_params[0] is not None: 
                    gamma_H2 = H2_params[0]
                else: 
                    gamma_H2 = DEFAULT_GAMMA
                if H2_params[1] is not None: 
                    n_H2 = H2_params[1]
                else: 
                    n_H2 = DEFAULT_N
            else: #if cannot find params by J, set them to default
                gamma_H2 = DEFAULT_GAMMA
                n_H2 = DEFAULT_N
        
            #get He params
            if He_dict.get(str(J_lower) + '_' + str(K_lower)) is not None: 
                He_params = He_dict.get(str(J_lower) + '_' + str(K_lower))
                if He_params[0] is not None: 
                    gamma_He = He_params[0]
                else: 
                    gamma_He = DEFAULT_GAMMA
                if He_params[1] is not None: 
                    n_He = He_params[1]
                else: 
                    n_He = DEFAULT_N
            elif He_dict.get(str(J_lower)) is not None: 
                He_params = He_dict.get(str(J_lower))
                if He_params[0] is not None: 
                    gamma_He = He_params[0]
                else: 
                    gamma_He = DEFAULT_GAMMA
                if He_params[1] is not None: 
                    n_He = He_params[1]
                else: 
                    n_He = DEFAULT_N
            else: #if cannot find params by J, set them to default
                gamma_He = DEFAULT_GAMMA
                n_He = DEFAULT_N
        
        
        v_ij = E_upper - E_lower
        
        data = [v_ij, A, E_lower, g_upper, gamma_H2, n_H2, gamma_He, n_He]
            
        for item in data: 
            f.write("%s " % item)
        f.write("\n")
        
        global counter
        counter += 1
        #print(counter)
    
    f.close()
    
    print("Write infile in %s seconds" % (time.time() - file_time))

    load_time = time.time()
    
    print("Bulk inserting exomol data...")

    cursor.execute("LOAD DATA LOCAL INFILE '/home/toma/Desktop/exomol.txt' INTO TABLE transitions FIELDS TERMINATED BY ' ' LINES TERMINATED BY '\n' \
              (@col1, @col2, @col3, @col4, @col5, @col6, @col7, @col8) SET nu=@col1, A=@col2, elower=@col3, g_upper=@col4, \
              gamma_H2=@col5, n_H2=@col6, gamma_He=@col7, n_He=@col8, line_source_id={}, particle_id={};".format(line_source_id, particle_id))
    
    print('Executed {} lines of exomol data'.format(counter))
    print("Bulk inserted exomol data in %s seconds" % (time.time() - load_time))    

################

start_time = time.time()

#connect to the database
db = MySQLdb.connect(host='localhost', user='toma', passwd='Happy810@', db='linelist') 

#create a cursor object
cursor = db.cursor()

#disable autocommit to improve performance

sql_order('SET autocommit = 0')
sql_order('SET unique_checks = 0')
sql_order('SET foreign_key_checks = 0')
sql_order('SET sql_log_bin = 0')

####################
    
#for CO
#insert_partitions('/home/toma/Desktop/12C-16O__Li2015_partition.pf', particle_id)
insert_partitions('/home/toma/Desktop/linelists-database/PH3_partitions.txt', 'EXOMOL_Li2015', particle_id)

###################

states_time = time.time()
#get parameters needed to insert exomol data into transitions
print('Loading huge ass states file')
#states in id order starts in 1 
#Es, gs, Js = np.loadtxt('/home/toma/Desktop/12C-16O__Li2015.states', usecols=(1, 2, 3), unpack=True)
Es, gs, Js, Ks= np.loadtxt('/home/toma/Desktop/linelists-database/PH3_states.txt', usecols=(1, 2, 3, 6), unpack=True)
print('Finished loading states file in %s seconds' % (time.time() - states_time))

H2_dict = temp_broad_param_dict('/home/toma/Desktop/linelists-database/PH3_H2_broad.txt')

He_dict = temp_broad_param_dict('/home/toma/Desktop/linelists-database/PH3_He_broad.txt')
'''
H2_dict = temp_broad_param_dict('/home/toma/Desktop/12C-16O__H2.broad')
He_dict = temp_broad_param_dict('/home/toma/Desktop/12C-16O__He.broad')
'''

#transition file
counter = 0 

#length_trans = sum(1 for line in open('/home/toma/Desktop/12C-16O__Li2015.trans'))

#with open('/home/toma/Desktop/12C-16O__Li2015.trans') as trans: #for CO
for file_num in range(1, 101):
    curr_file = '/home/toma/Desktop/linelists-database/PH3_trans_{}'.format(file_num)
    #curr_file = '/home/toma/Desktop/12C-16O__Li2015.trans'
    
    #get the number of lines in trans file
    length_trans = sum(1 for line in open(curr_file))
    print(length_trans, 'lines : Opened the transition file')
    
    #sql_order('ALTER TABLE transitions AUTO_INCREMENT = 1') #this can be used to solve auto_increment problems as well
    
    with open(curr_file) as trans:
        #for spliiting file into smalller chunks...but mysql auto_increment seems to not be working properly
        start_line = 0
        max_size = 1e7
        
        repeat = 0
        while length_trans >= start_line + max_size: 
            
            #need to change exomol_Li2015 into line_source_id next time running !!!!!!!!!!!!!
            insert_exomol(start_line, int(start_line + max_size), '/home/toma/Desktop/exomol.txt', trans, 'EXOMOL_Li2015', particle_id)
            
            #islice removes starts from the next line after the last read line
            length_trans -= max_size
            #print(int(length_trans))
            repeat += 1
        
        #out of the while loop when difference between start_line and the max lines in trans file is less than 1e6
        insert_exomol(start_line, int(length_trans), '/home/toma/Desktop/exomol.txt', trans, 'EXOMOL_Li2015', particle_id)
        
        #insert_exomol('/home/toma/Desktop/exomol.txt', trans, particle_id)
        
    #commit one file altogether at one time    
    db.commit()
    trans.close()
    
    print('Finished loading {} with {} lines of data'.format(curr_file, int(length_trans + repeat * max_size)))
    #print('Finished loading {} with {} lines of data'.format(curr_file, int(length_trans)))
    
    
    #set @id:=0; update mytable set id = (@id := @id + 1) order by id; for correcting auto_increment if needed
    
###################

#turn them back on

sql_order('SET unique_checks = 1')
sql_order('SET foreign_key_checks = 1')
sql_order('SET sql_log_bin = 1')

cursor.close()
db.close()

print("Finished in %s seconds" % (time.time() - start_time))

#################################
#################################
#################################
#final boss function

#fp stands for filepath
def import_exomol_data(mol_name, iso_name, version_name, trans_fp, states_fp, partitions_fp, \
                       broad_H2_fp=None, broad_He_fp=None, DEFAULT_GAMMA, DEFAULT_N, trans_file_num, reference_link):
    
    ###################
    
    global DEFAULT_GAMMA, DEFAULT_N
    DEFAULT_GAMMA = DEFAULT_GAMMA
    DEFAULT_N = DEFAULT_N   
    if DEFAULT_GAMMA is None or DEFAULT_N is None: 
        raise Exception('DEFAULT_GAMMA or DEFAULT_N should not be null')
        
    ###################
    
    one_iso_time = time.time()  
    #connect to the database
    db = MySQLdb.connect(host='localhost', user='toma', passwd='Happy810@', db='linelist') 
    #create a cursor object
    cursor = db.cursor()
    #disable autocommit to improve performance
    sql_order('SET autocommit = 0')
    sql_order('SET unique_checks = 0')
    sql_order('SET foreign_key_checks = 0')
    sql_order('SET sql_log_bin = 0')
    
    ##################
    
    #load H2/He params and in the mean while
    #insert the line_source into source_properties and get line_source_id
    if broad_H2_fp is None and broad_He_fp is None: #when no .broad files in exomol
        no_broadening_param = True
        
        insert_version_query = "INSERT INTO source_properties(line_source, max_temperature, max_nu, num_lines, bool_air, \
        bool_H2, bool_He, reference_link, line_source_id) VALUES('%s', null, null, null, 'NO', 'NO', 'NO', '%s', null);" % \
        (version_name, reference_link)
        
    elif broad_H2_fp is not None and broad_He_fp is not None: #when both H2 and He .broad files in exomol
        no_broadening_param = False
        H2_dict = temp_broad_param_dict(broad_H2_fp)
        He_dict = temp_broad_param_dict(broad_He_fp)
        
        insert_version_query = "INSERT INTO source_properties(line_source, max_temperature, max_nu, num_lines, bool_air, \
        bool_H2, bool_He, reference_link, line_source_id) VALUES('%s', null, null, null, 'NO', 'YES', 'YES', '%s', null);" % \
        (version_name, reference_link)
        
    else: 
        raise Exception('Should have either neither or both of the H2 and He broad param files')
        
    sql_order(insert_version_query)
    
    get_line_source_id_query = "SELECT line_source_id FROM source_properties WHERE line_source = '{}'".format(version_name)
    
    data = fetch(get_line_source_id_query)

    if len(data) != 1:
        raise Exception('should have exactly one line_source_id corresponding to one line_source')
        
    line_source_id = data[0][0]
    
    #####################
    
    #get particle_id
    get_particle_id = "SELECT particle_id FROM particles WHERE iso_name = '{}'".format(iso_name)
    data = fetch(get_particle_id)
    if len(data) != 1:
        raise Exception('iso_name should correspond to exactly one isotopologue in the database')
    particle_id = data[0][0]
    
    #insert partitions
    insert_partitions(partitions_fp, line_source_id, particle_id)
    
    #load states
    states_time = time.time()
    #get parameters needed to insert exomol data into transitions
    print('Loading huge ass states file')
    #states in id order starts in 1 
    Es, gs, Js, Ks= np.loadtxt(states_fp, usecols=(1, 2, 3, 6), unpack=True)
    print('Finished loading states file in %s seconds' % (time.time() - states_time))
    
    ######################
    
    #insert transition files
    global counter
    counter = 0 
    
    for file_num in range(1, trans_file_num + 1):
        curr_file = '/home/toma/Desktop/linelists-database/' + trans_fp + str(file_num)        
        #get the number of lines in trans file
        length_trans = sum(1 for line in open(curr_file))
        print(length_trans, 'lines : Opened the transition file')
        
        #sql_order('ALTER TABLE transitions AUTO_INCREMENT = 1') #this can be used to solve auto_increment problems as well
        
        with open(curr_file) as trans:
            #for spliiting file into smalller chunks...but mysql auto_increment seems to not be working properly
            start_line = 0
            max_size = 1e7            
            repeat = 0
            
            while length_trans >= start_line + max_size: 
                insert_exomol(start_line, int(start_line + max_size), '/home/toma/Desktop/exomol.txt', trans, line_source_id, particle_id, no_broadening_param)
                #islice starts from the next line after the last read line
                length_trans -= max_size
                #print(int(length_trans))
                repeat += 1
            
            #out of the while loop when difference between start_line and the max lines in trans file is less than max_size
            insert_exomol(start_line, int(length_trans), '/home/toma/Desktop/exomol.txt', trans, line_source_id, particle_id, no_broadening_param)
            
        #commit one file altogether at one time
        db.commit()
        trans.close()
        
        print('Finished loading {} with {} lines of data'.format(curr_file, int(length_trans + repeat * max_size)))        
        
        #set @id:=0; update mytable set id = (@id := @id + 1) order by id; for correcting auto_increment if needed
    
    #turn them back on
    sql_order('SET unique_checks = 1')
    sql_order('SET foreign_key_checks = 1')
    sql_order('SET sql_log_bin = 1')
    
    cursor.close()
    db.close()
    
    print("Finished inserting", counter, "lines of exomol", version_name, "data for", iso_name, "in %s seconds" % (time.time() - one_iso_time))
    