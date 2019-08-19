#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 13:16:47 2019

@author: toma
"""

import compute_absorption

#Interactive User Interface for Absorption Cross Section Calculator
#Working but need to refine more to be more complete and userr friendly

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
        
        output = compute_absorption.new_compute_all(v, T, p, isotopologue, line_source)
        print('The absorption cross section for {} in {} at {}, {}, {} is {} \n'.format(isotopologue, line_source, v, T, p, output))
        toma = input('Do you want to obtain another absorption cross section data, Sir? \nEnter: Y or n \n')
    elif toma.lower() == 'n':
        print('Farewell, Sir. May the flame guide you~ \n')
        break
    else: 
        toma = input('Please enter Y or n \n')
quit()