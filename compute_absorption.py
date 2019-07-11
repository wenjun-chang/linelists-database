#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 11:14:43 2019

@author: toma
"""

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
from query_functions import fetch
import numpy as np
import time
from astropy.modeling.models import Voigt1D
import numexpr as ne
######################

#store the value of speed of light in cm/s
c = 2.99792458e10
#planck constant erg*s
h = 6.62606885e-27
#boltzmann constant
k_B = 1.38064852e-16
#reference Temperature in K
T_ref = 296 
#store the value of c_2 = h * c / k_B
c_2 = h * c / k_B
#store the conversion from gram to amu
G_TO_AMU = 1.66053904e-24

########################

#given input v(nu), T, p, iso, source, and version

#fetch the partition function value given an input T, temperature
def get_partition(T): #temp has to be a float i.g. 19.0
    
    #query for the partition function given T, temperature
    query = "SELECT `partition` FROM partitions WHERE temperature = {}".format(T)
    
    data = fetch(query)
    
    if len(data) != 1:
        raise Exception('should have exactly one partition value given a specific T')
        
    return data[0][0]

#######################

#get the particle id and the isotopologue abundance of the input iso_name 
#iso_name format for example, CO2, is (13C)(16O)2
def get_particle(iso_name):
    query = "SELECT particle_id, iso_abundance, iso_mass FROM particles WHERE iso_name = '{}'".format(iso_name)
    
    data = fetch(query)
    
    if len(data) != 1:
        raise Exception('should have exactly one row for a specific particle')
    
    #data[0] = (particle_id, iso_abundance, iso_mass)
    return data[0]
    
##########################

#absorption(v, T, p) = S_ij(T) * f(v, v_ij, T, p)        
def compute_one_absorption(line, v, T, p, Q, iso_abundance, iso_mass):
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
    
    #compute line intensity function S_ij(T)
    S_ij = (iso_abundance * a * gp * math.exp(-c_2 * elower / T) * (1 - math.exp(-c_2 * v_ij / T))) / (8 * math.pi * c * math.pow(v_ij, 2) * Q)
    
    #compute gamma(p,T) for f
    #T_red = 296 K
    #gamma_p_T = p * ((T_ref / T)^n_H2 * gamma_H2 * f_H2 + (T_ref / T)^n_He * gamma_He * f_He)
    #where f_H2 = 0.85 and f_He = 0.15
    #if either n_H2 or n_He does not exist, f_H2/He (the exisiting one) = 1.0
    if n_H2 is not None and gamma_H2 is not None and n_He is not None and gamma_He is not None:
        gamma_p_T = p * (math.pow(T_ref/ T, n_H2) * gamma_H2 * 0.85 + math.pow(T_ref / T,n_He) * gamma_He * 0.15)
        
    #if n_H2 does not exist, f_He = 1
    elif (n_H2 is None  or gamma_H2 is None) and (n_He is not None and gamma_He is not None):
        gamma_p_T = p * math.pow(T_ref / T, n_He) * gamma_He
    
    #if n_He does not exist, f_H2 = 1
    elif (n_He is None or gamma_He is None) and (n_H2 is not None and gamma_H2 is not None):
        gamma_p_T = p * math.pow(T_ref / T, n_H2) * gamma_H2
    
    #if both n_H2 or n_He does not exist
    #gamma_p_T = p * (T_ref / T)^n_air * gamma_air
    else:
        gamma_p_T = p * math.pow(T_ref / T, n_air) * gamma_air
    

    #compute v_ij_star for f 
    #v_ij_star = v_ij + delta_net * p, where delta_net is computed in similar fashion to gamma_p_T
    if delta_H2 is not None and delta_He is not None: #if both delta exists
        v_ij_star = v_ij + p * (delta_H2 * 0.85 + delta_He * 0.15)
        
    elif delta_H2 is None and delta_He is not None: #when delta_H2 does not exist, f_He = 1.0
        v_ij_star = v_ij + p * delta_He
    
    elif delta_He is None and delta_H2 is not None: #when delta_He does not exist, f_H2 = 1.0
        v_ij_star = v_ij + p * delta_H2
        
    elif delta_air is not None: #when both delta_H2 and delta_He does not exist, use delta_air
        v_ij_star = v_ij + p * delta_air
    else: #when all deltas do not exist
        v_ij_star = v_ij
        
        
    #compute normalized line shape function f(v, v_ij, T, p)
    #############old 
    '''f = gamma_p_T / (math.pi * (math.pow(gamma_p_T, 2) + (v - v_ij_star)**2))'''
    
    
    
    #what is f_0
    droppler_broad = math.sqrt((8 * k_B * T * math.log(2)) / (iso_mass * G_TO_AMU * c**2)) * v_ij
    
    absorption_function = Voigt1D(x_0=v_ij_star, amplitude_L= 1 / (gamma_p_T * math.pi), fwhm_L=2 * gamma_p_T, fwhm_G=droppler_broad)
    
    absorption = S_ij * absorption_function(v)
    
    '''
    #compute absorption cross section
    absorption = S_ij * f
    '''
    return absorption

###################
    
#given input v, T, p, iso_name, source, and version, fetch all the line data of the input iso_name
#use the parameters feteched to compute absorption cross section with the help of other functions
def compute_all(v, T, p, iso_name, line_source, default=False): 
    
    #get paritition using the correct function
    Q = get_partition(T)
    
    #get particle_id and iso_abundance using the correct function
    particle_data = get_particle(iso_name)
    particle_id = particle_data[0]
    iso_abundance = particle_data[1]
    iso_mass = particle_data[2]

     #connect to the database
    db = MySQLdb.connect(host='localhost', user='toma', passwd='Happy810@', db='linelist') 
    #do put actual password when run
    
    #create a cursor object
    cursor = db.cursor()
    
    
    #####################add when default
    
    if default == True: 
        #use this query for getting the lines of default line source
        query = "SELECT nu, A, gamma_air, n_air, delta_air, elower, gp, gamma_H2, n_H2, delta_H2, gamma_He, n_He, \
        delta_He FROM transitions INNER JOIN particles ON particles.default_line_source = transitions.line_source \
        WHERE particles.particle_id = {};".format(particle_id)
    
    else:
        #query for all the lines of the specified isotopologue from the user given nu, line_source
        query = "SELECT nu, A, gamma_air, n_air, delta_air, elower, gp, gamma_H2, \
        n_H2, delta_H2, gamma_He, n_He, delta_He FROM transitions WHERE particle_id = {} AND \
        line_source = '{}'".format(particle_id, line_source)
    
    print(query)
    
    #this gives us a table of all the parameters we desire in a table in mysql
    cursor.execute(query)
    #rowcount is a read-only attribute and returns the number of rows that were affected by the execute() method.
    rows = cursor.rowcount
    print(rows)
    #the table could be gigantic, therefore fetchall() could be slow, therefore
    #would rather fetch one single line as a tuple ( , , , ) each time and
    #compute the absorption for that line, store it in variable and sum it over iterations. 
    absorption_cross_section = np.zeros(len(v))
    for i in range(rows):
        #fetch one line
        line = cursor.fetchone()
        cond =  np.logical_and(v >= line[0] - 25, v <= line[0] + 25)
        if np.sum(cond) > 0:
            absorption_cross_section[cond] += compute_one_absorption(line, v[cond], T, p, Q, iso_abundance, iso_mass)
        print(i)
        '''
        #use the line and given input to compute absorption and put it into a list
        for j in range(len(v)):
            if line[0] >= v[j] - 25 and line[0] <= v[j] + 25:             
                absorption = compute_one_absorption(line, v[j], T, p, Q, iso_abundance, iso_mass)
                absorption_cross_section[j] += absorption
        '''
    #close up cursor and connection
    cursor.close()
    db.close()
    
    return absorption_cross_section
        
#########################

#obtain input v, T, p, iso_name, line_source from the user
def main():
    #1,2,3,4,5
    #1,10,100,1000,10000
    #1,2,4,8,16

    start_time = time.time()
    
    #wavelengths = np.logspace(np.log10(0.3e-4), np.log10(30e-4), 4616)
    wavelengths = np.load('/home/toma/Desktop/wavelengths.npy') * 100
    wavenums = 1.0/wavelengths
    #wavelengths = np.exp(np.linspace(np.log(0.3e-4), np.log(30e-4), 4616))
    
    absorption_cross_section = compute_all(wavenums, 1000, 0.1, '(12C)(16O)', 'EXOMOL_Li2015')
    #print('absorption_cross_section is', absorption_cross_section)
    np.save('/home/toma/Desktop/absorption.npy', absorption_cross_section)
    
    print("Finished in %s seconds" % (time.time() - start_time))
    
    '''
    print('Hello! Welcome to Toma\'s linelist database, sir.')
    toma = input('Do you want to obtain absorption cross section data today, Sir? \nEnter: Y or n \n')
    while True: 
        if toma.lower() == 'y':
            isotopologue = input('What isotopologue of a molecule are you looking for? Format example: CO2 = (12C)(16O)2 \n')
            source = input('What type of data do Sir wish to compute from? \nEnter: HITRAN or EXOMOL or DEFAULT--the optimal'
                           'version for the molecule you chooose, Sir \n')
            if source.lower() != 'default': 
                hitran_versions = '2016, 2012, 2008'
                exomol_versions = 'Li2015, SAITY'
                if source.lower() == 'hitran':
                    version = input('What version of {} do Sir desire? \nVersions available include {}'
                                    ' (case sensitive) \n'.format(source.upper(), hitran_versions))
                if source.lower() == 'exomol':
                    version = input('What version of {} do Sir desire? \nVersions available include {}'
                                    ' (case sensitive) \n'.format(source.upper(), exomol_versions))
                    
                line_source = source.upper() + '_' + version
            else: #for default
                line_source = source.lower()
                
            wavenums = input('At what wavenumber(s) (cm^-1) Sir? \n')
            temp = wavenums.split(',')
            v = []
            for i in temp:
                v.append(float(i))
            T = float(input('At what temperature (K) Sir? \n'))
            p = float(input('At what pressure (atm) Sir? \n'))
            
            if line_source == 'default':
                output = compute_all(v, T, p, isotopologue, line_source, True)
            else: 
                output = compute_all(v, T, p, isotopologue, line_source)
            print('The absorption cross section for {} in {} at {}, {}, {} is {} \n'.format(isotopologue, line_source, v, T, p, output))
            toma = input('Do you want to obtain another absorption cross section data, Sir? \nEnter: Y or n \n')
        elif toma.lower() == 'n':
            print('Farewell, Sir. May the flame guide you~ \n')
            break
        else: 
            toma = input('Please enter Y or n \n')
    quit()
    '''

if __name__ == '__main__':
    main()