#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 13:12:27 2019

@author: toma
"""

#this code fetch data from database and compute absorption cross section data 
#as the user desired with a user interface main()

import MySQLdb
import math

########################

#given input v(nu), T, p, molecule, data_type, and version

#fetch the partition function value given an input T, temperature
def get_partition(T): #temp has to be a float i.g. 19.0
    #connect to the database
    db = MySQLdb.connect(host='localhost', user='toma', passwd='Happy810@', db='linelist') 
    #do put actual password when run
    
    #create a cursor object
    cursor = db.cursor()
    
    #query for the partition function given T, temperature
    query = '''select `partition` from partitions where temp = {}'''.format(T)
    
    try: 
        #there should only be one single value, otherwise the database is wrong
        cursor.execute(query)
        data = cursor.fetchall()
        
        #check if the database is correct
        if (data[1] == True):
            raise Exception('should only have one partition value given a specific T')
        
        partition_value = data[0][0]
        return partition_value
    
    except: 
        #if errors occur
        db.rollback()
        
    finally: 
        #close up cursor and connection
        cursor.close()
        db.close()

#######################

#get the isotope abundance of the input molecule  
#molecule_name format for example, CO2, is (13C)(16O)2
def get_particle(mol_name):
     #connect to the database
    db = MySQLdb.connect(host='localhost', user='toma', passwd='Happy810@', db='linelist')
    #do put actual password when run
    
    #create a cursor object
    cursor = db.cursor()
    
    #query for the partition function given T, temperature
    query = '''select iso_abundance from particles where molecule_name = {}'''.format(mol_name)
    
    try: 
        cursor.execute(query)
        data = cursor.fetchone()
        
        particle_id = data[0]
        iso_abundance = data[1]
        return (particle_id, iso_abundance)
    
    except: 
        #if errors occur
        db.rollback()
        
    finally: 
        #close up cursor and connection
        cursor.close()
        db.close()

##########################

#absorption(v, T, p) = S_ij(T) * f(v, v_ij, T, p)        
def compute_one_absorption(line, v, T, p, Q, iso_abundance):
    #line is a tuple returned by fetchone() operation
    #parameters for fectone() data corresponding to tuple indexes: 
    #(nu, a, gamma_air, n_air, delta_air, elower, gp, gamma_H2, n_H2, delta_H2, gamma_He, n_He, delta_He)
    #(0,  1,     2,       3,       4,        5,     6,    7,       8,      9,      10,      11,     12   )
    
    #just name the variables to make things clear and easy to check
    v_ij = line[0] #v_ij = nu !!!!!!!
    a = line[1]
    gamma_air = line[2]
    n_air = line[3]
    delta_air = line[4]
    elower = line[5]
    gp = line[6]
    gamma_H2 = line[7]
    n_H2 = line[8]
    delta_H2 = line[9]
    gamma_He = line[10]
    n_He = line[11]
    delta_He = line[12]
    
    #store the value of speed of light in m/s
    c = 299792458 
    #planck constant
    h = 6.62607004 * math.pow(10, -34)
    #boltzmann constant
    k_B = 1.38064852 * math.pow(10, -23)
    #store the value of c_2 = h * c / k_B
    c_2 = h * c / k_B
    
    #compute line intensity function S_ij(T)
    S_ij = (iso_abundance * a * gp * math.exp(-c_2 * elower / T) * (1 - math.exp(-c_2 * v_ij / T))) / (8 * math.pi * c * math.pow(v_ij, 2) * Q)
    
    
    #compute gamma(p,T) for f
    #T_red = 296 K
    #gamma_p_T = p * ((T_ref / T)^n_H2 * gamma_H2 * f_H2 + (T_ref / T)^n_He * gamma_He * f_He)
    #where f_H2 = 0.85 and f_He = 0.15
    #if either n_H2 or n_He does not exist, f_H2/He (the exisiting one) = 1.0
    if n_H2 != 0.0 and gamma_H2 != 0.0 and n_He != 0.0 and gamma_He != 0.0:
        gamma_p_T = p * (math.pow((296 / T), n_H2) * gamma_H2 * 0.85 + math.pow((296 / T),n_He) * gamma_He * 0.15)
        
    #if n_H2 does not exist, f_He = 1
    elif (n_H2 == 0.0 or gamma_H2 == 0.0) and (n_He != 0.0 and gamma_He != 0.0):
        gamma_p_T = p * math.pow((296 / T), n_He) * gamma_He
    
    #if n_He does not exist, f_H2 = 1
    elif (n_He == 0.0 or gamma_He == 0.0) and (n_H2 != 0.0 and gamma_H2 != 0.0):
        gamma_p_T = p * math.pow((296 / T), n_H2) * gamma_H2
    
    #if both n_H2 or n_He does not exist
    #gamma_p_T = p * (T_ref / T)^n_air * gamma_air
    else:
        gamma_p_T = p * math.pow((296 / T), n_air) * gamma_air
    

    #compute v_ij_star for f 
    #v_ij_star = v_ij + delta_net * p, where delta_net is computed in similar fashion to gamma_p_T
    if delta_H2 != 0.0 and delta_He != 0.0: #if both delta exists
        v_ij_star = v_ij + p * (delta_H2 * 0.85 + delta_He * 0.15)
        
    elif delta_H2 == 0.0 and delta_He != 0.0: #when delta_H2 does not exist, f_He = 1.0
        v_ij_star = v_ij + p * delta_He
    
    elif delta_He == 0.0 and delta_H2 != 0.0: #when delta_He does not exist, f_H2 = 1.0
        v_ij_star = v_ij + p * delta_H2
        
    else: #when both delta_H2 and delta_He does not exist, use delta_air
        v_ij_star = v_ij + p * delta_air
        
        
    #compute normalized line shape function f(v, v_ij, T, p)
    f = gamma_p_T / (math.pi * (math.pow(gamma_p_T, 2) + math.pow((v - v_ij_star), 2)))
    
    #compute absorption cross section
    absorption = S_ij * f
    return absorption

###################
    
#given input v, T, p, molecule, data_type, and version, fetch all the line data of the input molecule
#use the parameters feteched to compute absorption cross section with the help of other functions
def compute_all(v, T, p, molecule, data_type, version): 
    
    #get paritition using the correct function
    Q = get_partition(T)
    
    #get particle_id and iso_abundance using the correct function
    particle_data = get_particle(molecule)
    particle_id = particle_data[0]
    iso_abundance = particle_data[1]

     #connect to the database
    db = MySQLdb.connect(host='localhost', user='toma', passwd='Happy810@', db='linelist') 
    #do put actual password when run
    
    #create a cursor object
    cursor = db.cursor()
    
    #query for all the lines of the specified molecule from the user given nu, data_type and version
    query = '''select nu, a, gamma_air, n_air, delta_air, elower, gp, gamma_H2, \
    n_H2, delta_H2, gamma_He, n_He, delta_He from `lines` where particle_id = {} AND \
    data_type = {} AND version = {}'''.format(particle_id, data_type, version)
    
    try: 
        #this gives us a table of all the parameters we desire in a table in mysql
        cursor.execute(query)
        #rowcount is a read-only attribute and returns the number of rows that were affected by the execute() method.
        rows = cursor.rowcount
        #the table could be gigantic, therefore fetchall() could be slow, therefore
        #would rather fetch one single line as a tuple ( , , , ) each time and
        #compute the absorption for that line, store it in a list (just to be safe, not necessary), and sum it. 
        absorb_lst = []
        for i in range(len(rows)):
            #fetch one line
            line = cursor.fetchone()
            #use the line and given input to compute absorption and put it into a list
            absorption = compute_one_absorption(line, v, T, p, Q, iso_abundance)
            absorb_lst[i] = absorption
        
        #sum up the absorption from all lines
        absorption_cross_section = sum(absorb_lst)
        return absorption_cross_section
    
    except: 
        #if errors occur
        db.rollback()
        
    finally: 
        #close up cursor and connection
        cursor.close()
        db.close()
        
#########################

#obtain input v, T, p, molecule, data_type, and version from the user
def main():
    print('Hello! Welcome to Toma\'s linelist database, sir.')
    toma = input('Do you want to obtain absorption cross section data today, Sir? \n Enter: Y or n')
    while toma.lower == 'y':
        molecule = input('What molecule/isotope are you looking for? Format example: CO2 = (12C)(16O)2')
        data_type = input('What type of data do Sir wish to compute from? \n Enter: HITRAN or EXOMOL')
        version = input('What version of {} do Sir desire?'.format(data_type))
        v = input('At what wavenumber (cm^-1) Sir?')
        T = input('At what temperature (K) Sir?')
        p = input('At what pressure (atm) Sir?')
        
        output = compute_all(v, T, p, molecule, data_type, version)
        print('The absorption cross section for {} in {}, {} at {}, {}, {} is {}'.format(molecule, data_type, version, v, T, p, output))
        toma = input('Do you want to obtain another absorption cross section data, Sir? \n Enter: Y or n')
    quit()

if __name__ == '__main__':
    main()