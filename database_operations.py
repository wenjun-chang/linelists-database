#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 13:53:21 2019

@author: toma
"""
#this file imports all the package python files and create functions to do different operations
#if runned as main this file is going to compute absorption and plot the result
##################

import matplotlib.pyplot as plt
import numpy as np
from create_tables import create_database_with_tables_and_indexes
from hitran_auto_download import import_hitran_and_populate_particles_table
from hitemp_import import populate_all_hitemp
from exomol_auto_download import populate_all_exomol
from compute_absorption import new_compute_all
from plot_absorption import plot_absorption

##################

#this fucntion creates the database with indexes on nu, A, and line_source_id, 
#populates the particle table with all molecules in HITRAN, and populate the 
#transitions table with HITRAN (and HITEMP, optional) lines
#default to import both HITRAN and HITEMP, if dont want HITEMP pass in a False param
def create_database_and_populate_hitran_and_hitemp(import_hitemp=True): 
    create_database_with_tables_and_indexes() #takes ~5 secs
    import_hitran_and_populate_particles_table() #takes ~10 mins
    if import_hitemp is True: 
        populate_all_hitemp() #takes from 60 mins to 100 mins

#populates exomol data that are not done, not mega_file_molecules, and not 'H2S' or 'CaO' because these two molecules have some problem
#change the molecule list `done` and `mega_file_molecules` in exomol_auto_download.py
def populate_exomol():
    populate_all_exomol()

#cython, numba, and numexpr optimizations cannot beat the python array operation optimization
#this fucntion uses the numpy arrray operation one
def compute_absorption(output_fileaname, wavenum_file, T, p, iso_name, line_source='default'): 
    wavenums = np.loadtxt(wavenum_file, usecols=0, unpack=True) # for CO I guess
    absorption_cross_section = new_compute_all(wavenums, T, p, iso_name, line_source)
    out_fp = '/home/toma/Desktop/{}.npy'.format(output_fileaname)
    np.save(out_fp, absorption_cross_section)
    return out_fp

#plots the absorption corss section that was computed
def plot_calculated_absorption_against_wavelength(wavenum_file, absorption_file):
    plot_absorption(wavenum_file, absorption_file)

if __name__ == '__main__': 
    #absorption_file = compute_absorption('poster_demo', '/home/toma/Desktop/output_1000_1d-1.xsec', 1000, 0.1, '(12C)(16O)', 'HITRAN_2016')
    #plot_calculated_absorption_against_wavelength('/home/toma/Desktop/output_1000_1d-1.xsec', absorption_file)
    
    wavenums, absorptions = np.loadtxt('/home/toma/Desktop/output_1000_1d-1.xsec', usecols=(0, 1), unpack=True)
    wavelengths = (1 / wavenums) * 1e4
    plt.xscale('log')
    plt.yscale('log')
    plt.plot(wavelengths, absorptions, label='Absorption Plot')
    plt.gca().legend(('database absorption',), fontsize='xx-large')
    plt.title("Absorption cross section plot", fontsize='xx-large')
    plt.xlabel('Wavelength (um)', fontsize='xx-large')
    plt.ylabel('Absorption Cross Section (cm^2)', fontsize='xx-large')
    plt.grid(True)
    plt.show()