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

#numexpr
import numexpr as ne
ne.set_vml_accuracy_mode('fast')
ne.set_vml_num_threads(1)
import os; os.environ['NUMEXPR_MAX_THREADS'] = '1'

#numba
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
    
    #get line source id
    if line_source == 'default': 
        line_source_id = fetch("SELECT default_line_source_id FROM particles WHERE particle_id = {}".format(particle_id))[0][0]
    else: 
        get_line_source_id_query = "SELECT line_source_id FROM source_properties WHERE line_source = '{}' and \
        particle_id = {}".format(line_source, particle_id)
        data = fetch(get_line_source_id_query)
        if len(data) != 1:
            raise Exception('should have exactly one line_source_id corresponding to one line_source and isotopologue')   
        line_source_id = data[0][0]
    print(line_source_id)
        
    #if computing using hitemp data, use hitran partitions, so get hitran line_source_id for partitions
    if 'HITEMP' in line_source: 
        get_hitran_source_id_query = "SELECT line_source, line_source_id FROM source_properties WHERE particle_id = {}".format(particle_id)
        sources = fetch(get_hitran_source_id_query)
        hitran_id = -1
        for source in sources: 
            if source[0].startswith('HITRAN'): 
                hitran_id = source[1]
        if hitran_id == -1:
            raise Exception('This isotopologue has hitemp but no hitran linelist which is weird')
        #use hitran id to get partitions for hitemp
        Q = get_partition(T, hitran_id, particle_id)
    
    else: #for other sources, use line source id to get partitions   
        #get paritition using the correct function
        Q = get_partition(T, line_source_id, particle_id)
    print(Q)

    #connect to the database
    db = MySQLdb.connect(host='localhost', user='toma', passwd='Happy810@', db='linelist') 
    #do put actual password when run
    
    #create a cursor object
    cursor = db.cursor()
    
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
def compute_one_wavenum(wavenumber, T, p, iso_abundance, iso_mass, Q, v_ij_star, a, elower, g_upper, gamma_p_T):#, lower_indexes, upper_indexes):
    
    '''
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
    voigt_profile = np.real(wofz) / sigma_thermal / np.sqrt(2*pi)
    absorption = S_ij * voigt_profile
    
    #spectrum: T = e^(-absorption*n*l)
    #n = p / (Kb * T) 0.5/(1.38064852e-23*270)
    #n = 1 / (1.38064852e-16 * 270)
    #n = p / (1.38064852e-16 * T)
    #l = 0.9e9
    #absorption = np.exp(-absorption*n*l)
    #return np.exp(-np.sum(absorption)*n*l)
    
    return np.sum(absorption)
    
#########################
#@profile
def new_compute_all(v, T, p, iso_name, line_source='default'): 
    #get particle_id and iso_abundance using the correct function
    print(len(v))
    particle_data = get_particle(iso_name)
    particle_id = particle_data[0]
    iso_abundance = particle_data[1]
    iso_mass = particle_data[2]
    print(particle_id, iso_abundance, iso_mass)
    '''
    ##delete later###############use this temporarily for testing
    get_line_source_id_query = "SELECT line_source_id FROM source_properties WHERE line_source = '{}' and \
        particle_id = {}".format(line_source, particle_id)
    data = fetch(get_line_source_id_query)
    if len(data) != 1:
        raise Exception('should have exactly one line_source_id corresponding to one line_source and isotopologue')   
    line_source_id = data[0][0]
    print(line_source_id)
    
    '''
    #get line source id
    if line_source == 'default': 
        line_source_id = fetch("SELECT default_line_source_id FROM particles WHERE particle_id = {}".format(particle_id))[0][0]
    else: 
        get_line_source_id_query = "SELECT line_source_id FROM source_properties WHERE line_source = '{}' and \
        particle_id = {}".format(line_source, particle_id)
        data = fetch(get_line_source_id_query)
        if len(data) != 1:
            raise Exception('should have exactly one line_source_id corresponding to one line_source and isotopologue')   
        line_source_id = data[0][0]
    print(line_source_id)
        
    #if computing using hitemp data, use hitran partitions, so get hitran line_source_id for partitions
    if 'HITEMP' in line_source: 
        get_hitran_source_id_query = "SELECT line_source, line_source_id FROM source_properties WHERE particle_id = {}".format(particle_id)
        sources = fetch(get_hitran_source_id_query)
        hitran_id = -1
        for source in sources: 
            if source[0].startswith('HITRAN'): 
                hitran_id = source[1]
        if hitran_id == -1:
            raise Exception('This isotopologue has hitemp but no hitran linelist which is weird')
        #use hitran id to get partitions for hitemp
        Q = get_partition(T, hitran_id, particle_id)
    
    else: #for other sources, use line source id to get partitions   
        #get paritition using the correct function
        Q = get_partition(T, line_source_id, particle_id)
    print(Q)
    '''
    #Q = 162879.38910000 #NO2
    Q = 152.18884000 #H2O at 270K
    print(Q)
    '''
    fetch_time =  time.time()
    
    #connect to the database
    db = MySQLdb.connect(host='localhost', user='toma', passwd='Happy810@', db='linelist') 
    #do put actual password when run
    
    #create a cursor object
    cursor = db.cursor()

    #query for all the lines of the specified isotopologue from the user given nu, line_sources
    query = "SELECT nu, A, gamma_air, n_air, delta_air, elower, g_upper, gamma_H2, \
    n_H2, delta_H2, gamma_He, n_He, delta_He FROM transitions WHERE particle_id = {} AND \
    line_source_id = '{}' ORDER BY nu".format(particle_id, line_source_id)
    
    #this gives us a table of all the parameters we desire in a table in mysql
    cursor.execute(query)
    print('Finished querying and fetching line list in %s seconds' % (time.time() - fetch_time))
    '''
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
    num_rows = int(5e6)
    
    counter = 0
    '''
    #start = 0 
    #all_lines_array = np.load('/home/toma/Desktop/linelists-database/(14N)(16O)2 .npy')
    all_lines_array = np.loadtxt('/home/toma/Desktop/(1H)2(16O)_database_foramt', usecols=(0,1,2,3,4,5,6), skiprows=1)    
    print(all_lines_array.shape)
    
    #sort the lines
    #argsort = np.argsort(all_lines_array[:,0])
    #all_lines_array = all_lines_array[argsort]
    #print(len(all_lines_array))
    '''
    absorption_cross_section = np.zeros(len(v))
    while True: 
        counter += 1
        print(counter)
        
        #end = min(start + num_rows, len(all_lines_array))
        #lines_array = all_lines_array[start:end]
        #print(len(lines_array))

        lines_table = cursor.fetchmany(size=num_rows) #########################
        lines_array = np.asarray(lines_table, dtype=np.float32) ################
        print(lines_array.shape)
        
        ###############construct gamma and n arrays
        #assuming cond returns an array of indexes that are within the range of cut off
        #parameters for fetchall() data corresponding to tuple indexes: 
        #(nu, a, gamma_air, n_air, delta_air, elower, g_upper, gamma_H2, n_H2, delta_H2, gamma_He, n_He, delta_He)
        #(0,  1,     2,       3,       4,        5,     6,    7,       8,      9,      10,      11,     12   )
        
        v_ij = lines_array[:,0] #v_ij = nu !!!!!!!
        print(v_ij)
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
        gamma_p_T = p * (T_ref / T)**(n_air) * gamma_air
        v_ij_star = v_ij + p * delta_air
        
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
        #v_ij_star = v_ij + delta_net * p, wcounter += 1
        #here delta_net is computed in similar fashion to gamma_p_T
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
        
        for i in range(len(v)): 
            if i % 100000 == 0: 
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
        '''
    #close up cursor and connection
    #cursor.close()
    #db.close()
    print(counter)

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
    wavenums, test = np.loadtxt('/home/toma/Desktop/output_1000_1d-1.xsec', usecols=(0,1), unpack=True)
    
    #wavenums = np.loadtxt('/home/toma/Desktop/N2O_1000_1d-1.xsec', usecols=0, unpack=True)
    #wavelengths, continuum = np.loadtxt('/home/toma/Desktop/Telluric Spectrum', usecols=(2, 5) , skiprows=1, unpack=True)
    #wavenums = 1e4 / wavelengths
    
    #earth_spec = new_compute_all(wavenums, 270, 0.5, '(1H)2(16O)', 'HITRAN_2016')
    #earth_spec *= continuum
    #np.save('/home/toma/Desktop/telluric_spectrum_model.npy', earth_spec)
    #absorption_cross_section = new_compute_all(wavenums, 1000, 0.1, '(14N)(16O)2', 'HITEMP_2019')
    #absorption_cross_section = compute_all(wavenums, 1000, 0.1, '(14N)(16O)2', 'HITEMP_2019')
    absorption_cross_section = new_compute_all(wavenums, 1000, 0.1, '(12C)(16O)', 'HITRAN_2016')
    print('absorption_cross_section is', absorption_cross_section)
    print(test)
    print(wavenums)
    print(absorption_cross_section.shape, test.shape)
    np.save('/home/toma/Desktop/absorption.npy', absorption_cross_section)
    
    print("Finished in %s seconds" % (time.time() - start_time))
if __name__ == '__main__':
    main()