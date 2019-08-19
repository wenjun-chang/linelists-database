#cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, linetrace=True

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 11:42:48 2019

@author: toma
"""

#this code fetch data from database and compute absorption cross section data 
#as the user desired with a user interface main()
print('new version')
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

#using cython
from libc.math cimport exp as c_exp, sqrt, M_PI
cdef extern from "vfastexp.h":
    double exp_approx "EXP" (double)
cdef extern from "Faddeeva.hh":
    complex Faddeeva "Faddeeva::w" (complex, double)

cimport cython #for profiling


#M_PI = np.pi
#complex.real working
#c_exp 68 secs vs. exp_approx 42 secs

######################

#store the value of speed of light in cm/s
cdef double c = 2.99792458e10
#planck constant erg*s
cdef double h = 6.62606885e-27
#boltzmann constant
cdef double k_B = 1.38064852e-16
#reference Temperature in K
cdef double T_ref = 296 
#store the value of c_2 = h * c / k_B
cdef double c_2 = h * c / k_B
#store the conversion from gram to amu
cdef double G_TO_AMU = 1.66054e-24#1.66053904e-24
#store pi value
cdef double pi = 3.1415926535897932384626433

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
    
#########################

@cython.binding(True)
#@profile
def compute_one_wavenum(double wavenumber, double T, double p, double iso_abundance, double iso_mass, \
                        double Q, double[:] v_ij_star, double[:] a, double[:] elower, double[:] g_upper, \
                        double[:] gamma_p_T):
    cdef double S_ij, sigma_thermal, voigt_profile, absorption = 0
    cdef double complex z, wofz
    cdef int N = a.shape[0]
    cdef unsigned int i
    iso_mass = 46.0055
    iso_abundance = 1

    for i in range(N):
        S_ij = iso_abundance * a[i] * g_upper[i] * c_exp(-c_2 * elower[i] / T) * (1 - c_exp(-c_2 * v_ij_star[i] / T)) / (8 * M_PI * c * v_ij_star[i] * v_ij_star[i] * Q)
        sigma_thermal = sqrt(k_B * T / (iso_mass * G_TO_AMU * c * c)) * v_ij_star[i]
        z = (wavenumber - v_ij_star[i] + gamma_p_T[i] * 1j) / (sigma_thermal * sqrt(2))
        wofz = Faddeeva(z, 1)
        #wofz = scipy.special.wofz(z)
        voigt_profile = wofz.real / (sigma_thermal * sqrt(2 * M_PI))
        absorption += S_ij * voigt_profile
        
    return absorption

#########################
#@profile
@cython.binding(True)
cdef double[:] new_compute_all(double[:] v, double T, double p, iso_name, line_source='default'): 
    #get particle_id and iso_abundance using the correct function
    particle_data = get_particle(iso_name)
    cdef int particle_id
    cdef double iso_abundance, iso_mass
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
    '''
    cdef double Q = 162879.38910000 #NO2
    print(Q)
    '''
    fetch_time =  time.time()
    
    #connect to the database
    db = MySQLdb.connect(host='localhost', user='toma', passwd='Happy810@', db='linelist') 
    #do put actual password when run
    
    #create a cursor object
    cursor = db.cursor()
    
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
    cdef unsigned int num_rows = 5000000 #5e6
    cdef unsigned int start = 0, counter = 0
    cdef unsigned int length_lines_array,  N = len(v), i
    cdef double[:] absorption_cross_section = np.zeros(len(v))

    all_lines_array = np.array(np.load('/home/toma/Desktop/linelists-database/(14N)(16O)2 .npy'), dtype=np.float64)
    print(len(all_lines_array))
    
    while True: 
        counter += 1
        print(counter)
        
        end = min(start + num_rows, len(all_lines_array))
        lines_array = all_lines_array[start:end]
        
        length_lines_array = len(lines_array)
        print(length_lines_array)

        #lines_table = cursor.fetchmany(size=num_rows) #########################
        #lines_array = np.asarray(lines_table, dtype=np.float32) ################3
        
        ###############construct gamma and n arrays
        #assuming cond returns an array of indexes that are within the range of cut off
        #parameters for fetchall() data corresponding to tuple indexes: 
        #(nu, a, gamma_air, n_air, delta_air, elower, g_upper, gamma_H2, n_H2, delta_H2, gamma_He, n_He, delta_He)
        #(0,  1,     2,       3,       4,        5,     6,    7,       8,      9,      10,      11,     12   )
        '''
        size = len(lines_array)
        cdef double[:] v_ij[size], a[size], gamma_air[size], n_air[size], delta_air[size], elower[size], gamma_H2[size], \
        n_H2[size], delta_H2[size], gamma_He[size], n_He[size], delta_He[size]
        cdef int g_upper[size]
        '''
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
        assert(np.all(np.sort(lines_array[:,0]) == lines_array[:,0]))

        lower_indexes = np.searchsorted(lines_array[:,0], np.array(v) - 25, side='right') #where lines_array[indexes - 1] <= v - 25 < lines_array[indexes]
        upper_indexes = np.searchsorted(lines_array[:,0], np.array(v) + 25) #where lines_array[indexes - 1] < v + 25 <= lines_array[indexes]
        print(lower_indexes, upper_indexes)
        
        for i in range(N):
        #for i in range(len(v)): 
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
            if i == 0:
                print(v_ij_star[0], a[0], \
                                    elower[0], g_upper[0], \
                                    gamma_p_T[0])
        if length_lines_array < num_rows:
            break

    #close up cursor and connection
    #cursor.close()
    #db.close()
    
    return absorption_cross_section
    
#########################

#obtain input v, T, p, iso_name, line_source from the user
@cython.binding(True)
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
#if __name__ == '__main__':
#    main()
