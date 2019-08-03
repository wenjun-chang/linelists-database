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
import scipy

import numexpr as ne
ne.set_vml_accuracy_mode('fast')
ne.set_vml_num_threads(1)
import os; os.environ['NUMEXPR_MAX_THREADS'] = '1'

from numba import jit


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
G_TO_AMU = 1.66054e-24#1.66053904e-24
#store pi value
pi = 3.1415926535897932384626433

########################

#given input v(nu), T, p, iso, source, and version

#fetch the partition function value given an input T, temperature
def get_partition(T, line_source_id, particle_id): #temp has to be a float i.g. 19.0
    print(T, line_source_id, particle_id)
    #query for the partition function given T, temperature
    query = "SELECT `partition` FROM partitions WHERE temperature = {} AND line_source_id = {} AND particle_id = {}".format(T, line_source_id, particle_id)
    
    data = fetch(query)
    
    if len(data) != 1:
        raise Exception('should have exactly one partition value given a specific T and line source')
        
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
    #(nu, a, gamma_air, n_air, delta_air, elower, g_upper, gamma_H2, n_H2, delta_H2, gamma_He, n_He, delta_He)
    #(0,  1,     2,       3,       4,        5,     6,    7,       8,      9,      10,      11,     12   )
    
    #just name the variables to make things clear and easy to check
    v_ij = line[0] #v_ij = nu !!!!!!!
    a = line[1]
    gamma_air = line[2]
    n_air = line[3]
    delta_air = line[4]
    elower = line[5]
    g_upper = line[6]
    gamma_H2 = line[7]
    n_H2 = line[8]
    delta_H2 = line[9]
    gamma_He = line[10]
    n_He = line[11]
    delta_He = line[12]
    
    #compute line intensity function S_ij(T)
    S_ij = (iso_abundance * a * g_upper * math.exp(-c_2 * elower / T) * (1 - math.exp(-c_2 * v_ij / T))) / (8 * math.pi * c * math.pow(v_ij, 2) * Q)
    
    #compute gamma(p,T) for f
    #T_red = 296 K
    #gamma_p_T = p * ((T_ref / T)^n_H2 * gamma_H2 * f_H2 + (T_ref / T)^n_He * gamma_He * f_He)
    #where f_H2 = 0.85 and f_He = 0.15
    #if either n_H2 or n_He does not exist, f_H2/He (the exisiting one) = 1.0
    if n_H2 is not None and gamma_H2 is not None and n_He is not None and gamma_He is not None:
        gamma_p_T = p * (math.pow(T_ref/ T, n_H2) * gamma_H2 * 0.85 + math.pow(T_ref / T, n_He) * gamma_He * 0.15)
        
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
    
    ######################!!!!!!!!!!!!!!this statement might be problematic..check logic!
    elif delta_air is not None: #when both delta_H2 and delta_He does not exist, use delta_air
        v_ij_star = v_ij + p * delta_air
    else: #when all deltas do not exist
        v_ij_star = v_ij
    
    #alternative way to compute voigt function: 70.0055468082428 seconds ; 0.01% diff ; 18% speed up
    sigma_thermal = np.sqrt(k_B * T / iso_mass / G_TO_AMU / c**2) * v_ij_star
    z = (v - v_ij_star + gamma_p_T * 1j) / sigma_thermal / np.sqrt(2)
    voigt_profile = np.real(scipy.special.wofz(z))/sigma_thermal / np.sqrt(2*math.pi)
    absorption = S_ij * voigt_profile
    
    '''
    #astropy way computing voigt function : 86.21910214424133 seconds ; 0.01% diff
    doppler_broad = math.sqrt((8 * k_B * T * math.log(2)) / (iso_mass * G_TO_AMU * c**2)) * v_ij
    absorption_function = Voigt1D(x_0=v_ij_star, amplitude_L= 1 / (gamma_p_T * math.pi), fwhm_L=2 * gamma_p_T, fwhm_G=doppler_broad)
    absorption = S_ij * absorption_function(v)
    '''
    
    return absorption

###################
    
#given input v, T, p, iso_name, source, and version, fetch all the line data of the input iso_name
#use the parameters feteched to compute absorption cross section with the help of other functions
def compute_all(v, T, p, iso_name, line_source='default'): 
    
    #get particle_id and iso_abundance using the correct function
    particle_data = get_particle(iso_name)
    particle_id = particle_data[0]
    iso_abundance = particle_data[1]
    iso_mass = particle_data[2]

    #get paritition using the correct function
    if 'HITRAN' not in line_source: 
        
        get_line_source_id_query = "SELECT line_source_id FROM source_properties WHERE line_source = '{}' and \
        particle_id = {}".format(line_source, particle_id)
        data = fetch(get_line_source_id_query)
        if len(data) != 1:
            raise Exception('should have exactly one line_source_id corresponding to one line_source and isotopologue')   
        line_source_id = data[0][0]
        
        Q = get_partition(T, line_source_id, particle_id)
    else: #if computing using hitran data
        raise Exception('WAAAAAAAAAAAAA need to create the dictionary for best partition fucntion first')
        #or just randomly choose a version
        #Q = fetch("SELECT `partition` FROM partitions WHERE temperature = {} AND particle_id = {}".format(T, particle_id))[0][0]
        
    #connect to the database
    db = MySQLdb.connect(host='localhost', user='toma', passwd='Happy810@', db='linelist') 
    #do put actual password when run
    
    #create a cursor object
    cursor = db.cursor()
    
    if line_source == 'default': 
        #use this query for getting the lines of default line source
        query = "SELECT nu, A, gamma_air, n_air, delta_air, elower, g_upper, gamma_H2, n_H2, delta_H2, gamma_He, n_He, \
        delta_He FROM transitions INNER JOIN particles ON particles.default_line_source_id = transitions.line_source_id \
        WHERE particles.particle_id = {};".format(particle_id)
    
    else:
        #query for all the lines of the specified isotopologue from the user given nu, line_source
        query = "SELECT nu, A, gamma_air, n_air, delta_air, elower, g_upper, gamma_H2, \
        n_H2, delta_H2, gamma_He, n_He, delta_He FROM transitions WHERE particle_id = {} AND \
        line_source_id = '{}'".format(particle_id, line_source_id)
    
    print(query)
    
    #this gives us a table of all the parameters we desire in a table in mysql
    cursor.execute(query)
    #rowcount is a read-only attribute and returns the number of rows that were affected by the execute() method.
    rows = cursor.rowcount
    print(rows, 'lines')
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
        #print(i)
    
    #close up cursor and connection
    cursor.close()
    db.close()
    
    return absorption_cross_section
        
#########################
#@@profile
#@jit(parallel=False, fastmath=True)
def compute_one_wavenum(wavenumber, T, p, iso_abundance, iso_mass, Q, v_ij_star, a, elower, g_upper, gamma_p_T, lower_indexes, upper_indexes):
    

    wavenumber = np.resize(wavenumber, (wavenumber.size, v_ij_star.size)).T
    for i in range(len(wavenumber)):
        wavenumber[i][:lower_indexes[i]] = np.nan
        wavenumber[i][upper_indexes[i] + 1:] = np.nan
    
    #alternative way to compute voigt function: 70.0055468082428 seconds ; 0.01% diff ; 18% speed up
    #sigma_thermal = np.sqrt(k_B * T / iso_mass / G_TO_AMU / c**2) * v_ij_star
    z = (wavenumber.T - v_ij_star + gamma_p_T * 1j) / (np.sqrt(k_B * T / iso_mass / G_TO_AMU / c**2) * v_ij_star) / np.sqrt(2)
    wofz = scipy.special.wofz(z)
    voigt_profile = np.real(wofz)/ (np.sqrt(k_B * T / iso_mass / G_TO_AMU / c**2) * v_ij_star) / np.sqrt(2*np.pi)
    
    bool_nan = voigt_profile != np.nan
    absorption = voigt_profile
    #compute line intensity function S_ij(T)
    #S_ij = (iso_abundance * a * g_upper * np.exp(-c_2 * elower / T) * (1 - np.exp(-c_2 * v_ij_star / T))) / (8 * np.pi * c * v_ij_star**2 * Q)
    absorption[bool_nan] *= (iso_abundance * a * g_upper * np.exp(-c_2 * elower / T) * (1 - np.exp(-c_2 * v_ij_star / T))) / (8 * np.pi * c * v_ij_star**2 * Q)
    absorption[~bool_nan] = 0
    
    return np.sum(absorption, axis=1)
    '''

    #compute line intensity function S_ij(T)
    S_ij = iso_abundance * a * g_upper * np.exp(-c_2 * elower / T) * (1 - np.exp(-c_2 * v_ij_star / T)) / (8 * pi * c * v_ij_star**2 * Q)
    
    #alternative way to compute voigt function: 70.0055468082428 seconds ; 0.01% diff ; 18% speed up
    sigma_thermal = np.sqrt(k_B * T / iso_mass / G_TO_AMU / c**2) * v_ij_star
    z = (wavenumber - v_ij_star + gamma_p_T * 1j) / sigma_thermal / np.sqrt(2)
    wofz = scipy.special.wofz(z)
    voigt_profile = np.real(wofz)/sigma_thermal / np.sqrt(2*pi)
    absorption = S_ij * voigt_profile
    
    return np.sum(absorption)
    '''
    
#########################
#@profile
def new_compute_all(v, T, p, iso_name, line_source='default'): 
    #get particle_id and iso_abundance using the correct function
    particle_data = get_particle(iso_name)
    particle_id = particle_data[0]
    iso_abundance = particle_data[1]
    iso_mass = particle_data[2]    
    
    #use this temporarily for testing
    get_line_source_id_query = "SELECT line_source_id FROM source_properties WHERE line_source = '{}' and \
        particle_id = {}".format(line_source, particle_id)
    data = fetch(get_line_source_id_query)
    if len(data) != 1:
        raise Exception('should have exactly one line_source_id corresponding to one line_source and isotopologue')   
    line_source_id = data[0][0]
    print(line_source_id)
    
    '''
    #get paritition using the correct function
    if 'HITRAN' not in line_source and 'HITEMP' not in line_source: 
        Q = get_partition(T, line_source_id, particle_id)
        
    else: #if computing using hitran data
        raise Exception('WAAAAAAAAAAAAA need to create the dictionary for best partition fucntion first')
        #or just randomly choose a version
        #Q = fetch("SELECT `partition` FROM partitions WHERE temperature = {} AND particle_id = {}".format(T, particle_id))[0][0]
    '''
    Q = 162879.38910000 #NO2
    print(Q)
    '''
    fetch_time =  time.time()
    
    #connect to the database
    db = MySQLdb.connect(host='localhost', user='toma', passwd='Happy810@', db='linelist') 
    #do put actual password when run
    
    #create a cursor object
    cursor = db.cursor()
    
    if line_source == 'default': 
        #use this query for getting the lines of default line source
        query = "SELECT nu, A, gamma_air, n_air, delta_air, elower, g_upper, gamma_H2, n_H2, delta_H2, gamma_He, n_He, \
        delta_He FROM transitions INNER JOIN particles ON particles.default_line_source_id = transitions.line_source_id \
        WHERE particles.particle_id = {};".format(particle_id)
    
    else:
        #query for all the lines of the3.7729922865792e-27 specified isotopologue from the user given nu, line_source
        query = "SELECT nu, A, gamma_air, n_air, delta_air, elower, g_upper, gamma_H2, \
        n_H2, delta_H2, gamma_He, n_He, delta_He FROM transitions WHERE particle_id = {} AND \
        line_source_id = '{}'".format(particle_id, line_source_id)
    
    #this gives us a table of all the parameters we desire in a table in mysql
    cursor.execute(query)
    print('Finished querying and fetching line list in %s seconds' % (time.time() - fetch_time))
    
    #ran out of ~10GB of memory
    lines_table = cursor.fetchall()
    lines_array = np.asarray(lines_table)
    
    print('Finished querying and fetching line list in %s seconds' % (time.time() - fetch_time))
    
    absorption_cross_section = np.zeros(len(v))
    
    print(len(v))
    for i in range(len(v)): 
        print(i)
        absorption_cross_section[i] = compute_one_wavenum(v[i], T, p, iso_abundance, iso_mass, Q, lines_array)
    '''
    num_rows = int(1e3)
    start = 0 
    counter = 0
    all_lines_array = np.load('/home/toma/Desktop/linelists-database/(14N)(16O)2 .npy')
    absorption_cross_section = np.zeros(len(v))
    while True: 
        counter += 1
        print(counter)
        
        end = min(start + num_rows, len(all_lines_array))
        lines_array = all_lines_array[start:end]
        print(len(all_lines_array))
        print(len(lines_array))

        #lines_table = cursor.fetchmany(size=num_rows) #########################
        #lines_array = np.asarray(lines_table, dtype=np.float32) ################3
        
        ###############construct gamma and n arrays
        #assuming cond returns an array of indexes that are within the range of cut off
        #parameters for fetchall() data corresponding to tuple indexes: 
        #(nu, a, gamma_air, n_air, delta_air, elower, g_upper, gamma_H2, n_H2, delta_H2, gamma_He, n_He, delta_He)
        #(0,  1,     2,       3,       4,        5,     6,    7,       8,      9,      10,      11,     12   )
        
        v_ij = lines_array[:,0] #v_ij = nu !!!!!!!
        a = lines_array[:,1]
        gamma_air = lines_array[:,2]
        n_air = lines_array[:,3]
        delta_air = lines_array[:,4]
        elower = lines_array[:,5]
        g_upper = lines_array[:,6]
        gamma_H2 = lines_array[:,7]
        n_H2 = lines_array[:,8]
        delta_H2 = lines_array[:,9]
        gamma_He = lines_array[:,10]
        n_He = lines_array[:,11]
        delta_He = lines_array[:,12]
        
        #################
        #initialize an np array for gamma_p_T
        gamma_p_T = np.zeros(len(a))
        
        #arrays that check whether that parameter is null in that index
        bool_gamma_H2 = ~np.isnan(gamma_H2)
        bool_gamma_He = ~np.isnan(gamma_He)
        bool_n_H2 = ~np.isnan(n_H2)
        bool_n_He = ~np.isnan(n_He)
        
        #compute gamma(p,T) for f
        #T_red = 296 K
        #gamma_p_T = p * ((T_ref / T)^n_H2 * gamma_H2 * f_H2 + (T_ref / T)^n_He * gamma_He * f_He)
        #where f_H2 = 0.85 and f_He = 0.15
        #if either n_H2 or n_He does not exist, f_H2/He (the exisiting one) = 1.0
        has_H2_and_He_gamma_N = np.all([bool_gamma_H2, bool_n_H2, bool_gamma_He, bool_n_He], axis=0)
        gamma_p_T[has_H2_and_He_gamma_N] = p * (T_ref/ T)**(n_H2[has_H2_and_He_gamma_N]) * gamma_H2[has_H2_and_He_gamma_N] \
                 * 0.85 + (T_ref / T)**(n_He[has_H2_and_He_gamma_N]) * gamma_He[has_H2_and_He_gamma_N] * 0.15
        
        #if n_H2 does not exist, f_He = 1
        has_He_but_not_H2_gamma_N = np.all([bool_gamma_He, bool_n_He, ~np.logical_or(bool_gamma_H2, bool_n_H2)], axis=0)
        gamma_p_T[has_He_but_not_H2_gamma_N] = p * (T_ref / T)**(n_He[has_He_but_not_H2_gamma_N]) * gamma_He[has_He_but_not_H2_gamma_N]
    
        #if n_He does not exist, f_H2 = 1
        has_H2_but_not_He_gamma_N = np.all([bool_gamma_H2, bool_n_H2, ~np.logical_or(bool_gamma_He, bool_n_He)], axis=0)
        gamma_p_T[has_H2_but_not_He_gamma_N] = p * (T_ref / T)**(n_H2[has_H2_but_not_He_gamma_N]) * gamma_H2[has_H2_but_not_He_gamma_N]
    
        #if both n_H2 or n_He does not exist
        #gamma_p_T = p * (T_ref / T)^n_air * gamma_air
        has_only_air_gamma_N = gamma_p_T == 0
        gamma_p_T[has_only_air_gamma_N] = p * (T_ref / T)**(n_air[has_only_air_gamma_N]) * gamma_air[has_only_air_gamma_N]
        
        ###################
        #initialize an np array for v_ij_star
        v_ij_star = np.zeros(len(a))
        
        #arrays that check whether that parameter is null in that index
        bool_delta_H2 = ~np.isnan(delta_H2)
        bool_delta_He = ~np.isnan(delta_He)
        bool_delta_air = ~np.isnan(delta_air)
         
        #compute v_ij_star for f 
        #v_ij_star = v_ij + delta_net * p, where delta_net is computed in similar fashion to gamma_p_T
        has_H2_and_He_delta = np.logical_and(bool_delta_H2, bool_delta_He)
        v_ij_star[has_H2_and_He_delta] = v_ij[has_H2_and_He_delta] + p * (delta_H2[has_H2_and_He_delta] * 0.85 + \
                 delta_He[has_H2_and_He_delta] * 0.15)
        
        #when delta_H2 does not exist, f_He = 1.0
        has_He_but_not_H2_delta = np.logical_and(~bool_delta_H2, bool_delta_He)
        v_ij_star[has_He_but_not_H2_delta] = v_ij[has_He_but_not_H2_delta] + p * delta_He[has_He_but_not_H2_delta]
        
        #when delta_He does not exist, f_H2 = 1.0
        has_H2_but_not_He_delta = np.logical_and(bool_delta_H2, ~bool_delta_He)
        v_ij_star[has_H2_but_not_He_delta] = v_ij[has_H2_but_not_He_delta] + p * delta_H2[has_H2_but_not_He_delta]
        
        #when both delta_H2 and delta_He does not exist, use delta_air
        has_air_but_not_H2_and_He_delta = np.all([bool_delta_air, ~bool_delta_H2, ~bool_delta_He], axis=0)
        v_ij_star[has_air_but_not_H2_and_He_delta] = v_ij[has_air_but_not_H2_and_He_delta] + p * delta_air[has_air_but_not_H2_and_He_delta]
        
        #when all deltas do not exist
        has_no_delta = np.all([~bool_delta_air, ~bool_delta_H2, ~bool_delta_He], axis=0)
        v_ij_star[has_no_delta] = v_ij[has_no_delta]
        
        #need to pass in: v_ij_star, a, elower, g_upper, gamma_p_T : all arrays
        
        ##################
        #indexes = np.searchsorted(ines_array, v) #where v[indexes - 1] < lines_array <= v[indexes]
        lower_indexes = np.searchsorted(lines_array[:,0], v - 25, side='right') #where lines_array[indexes - 1] <= v - 25 < lines_array[indexes]
        upper_indexes = np.searchsorted(lines_array[:,0], v + 25) #where lines_array[indexes - 1] < v + 25 <= lines_array[indexes]
        print(lower_indexes, upper_indexes)
        '''
        for i in range(len(v)): 
            if i % 10000 == 0: 
                print(i)
                
            #no_need_compute = lower_indexes == upper_indexes
            #print(no_need_compute)
            #lower_indexes = lower_indexes[no_need_compute]
            #print(lower_indexes)
            
            absorption_cross_section[i] = compute_one_wavenum(v[i], T, p, iso_abundance, iso_mass, Q, \
                                    v_ij_star[lower_indexes[i] : upper_indexes[i]], a[lower_indexes[i] : upper_indexes[i]], \
                                    elower[lower_indexes[i] : upper_indexes[i]], g_upper[lower_indexes[i] : upper_indexes[i]], \
                                    gamma_p_T[lower_indexes[i] : upper_indexes[i]])
        if len(lines_array) < num_rows:
            break
        '''
        absorption_cross_section += compute_one_wavenum(v, T, p, iso_abundance, iso_mass, Q, v_ij_star, a, elower, g_upper, gamma_p_T, lower_indexes, upper_indexes)
        start += num_rows
        if end == len(all_lines_array):
            break
        
    #close up cursor and connection
    #cursor.close()
    #db.close()
    
    return absorption_cross_section
    
#########################

#obtain input v, T, p, iso_name, line_source from the user
def main():
        
    start_time = time.time()
    
    #wavelengths in m --> convert to cm (x100)
    #wavelengths = np.loadtxt('/home/toma/Desktop/output_1000_1d-1.xsec', usecols=1, unpack=True) * 100
    #thus, wavenums in cm^-1
    #wavenums = 1.0/wavelengths
    #wavelengths = np.exp(np.linspace(np.log(0.3e-4), np.log(30e-4), 4616))
    #wavenums = np.loadtxt('/home/toma/Desktop/output_1000_1d-1.xsec', usecols=0, unpack=True)
    
    wavenums = np.loadtxt('/home/toma/Desktop/N2O_1000_1d-1.xsec', usecols=0, unpack=True)
    
    absorption_cross_section = new_compute_all(wavenums, 1000, 0.1, '(14N)(16O)2', 'HITEMP_2019')
    #absorption_cross_section = compute_all(wavenums, 1000, 0.1, '(14N)(16O)2', 'HITEMP_2019')
    #print('absorption_cross_section is', absorption_cross_section)
    np.save('/home/toma/Desktop/absorption.npy', absorption_cross_section)
    
    print("Finished in %s seconds" % (time.time() - start_time))
    
    '''
    #interactive user interface
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