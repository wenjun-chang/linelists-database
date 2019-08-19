#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 17:46:30 2019

@author: toma
"""

#insert correct default line source id instead of '1' when initially inserting all the molecules 
#exomol > hitemp > hitran

import time
from query_functions import fetch, sql_order

def correct_default_line_source_id():
    
    start_time = time.time()
    
    total_particle_id = 242
    for i in range(total_particle_id):
        particle_id = i + 1
        print(particle_id)
        #get all the sources for that isotopologue
        default_line_source_id = -1
        get_line_sources = "SELECT line_source, line_source_id FROM source_properties WHERE particle_id = {};".format(particle_id)
        sources = fetch(get_line_sources)
        print(sources)
        has_exomol = False
        exomol_id = -1
        has_hitemp = False
        hitemp_id = -1
        has_hitran = False
        hitran_id = -1
        for one_source in sources: 
            source_name = one_source[0]
            if source_name.startswith('EXOMOL'):
                has_exomol = True
                exomol_id = one_source[1]
            elif source_name.startswith('HITEMP'):
                has_hitemp = True
                hitemp_id = one_source[1]
            elif source_name.startswith('HITRAN'):
                has_hitran = True
                hitran_id = one_source[1]
        if not has_exomol and not has_hitemp and not has_hitran:
            raise Exception('Oh Damn this isotopologue has none of the versions HITRAN, HITEMP, or EXOMOL...umm problematic~')
        if has_exomol: 
            default_line_source_id = exomol_id
        elif not has_exomol and has_hitemp: 
            default_line_source_id = hitemp_id
        else: #only HITRAN
            default_line_source_id = hitran_id
        print(default_line_source_id)
        update_default_line_source_id = 'UPDATE particles SET default_line_source_id = {} WHERE particle_id = {};'.format(default_line_source_id, particle_id)
        sql_order(update_default_line_source_id)
        print('Finished correcting particle ' + str(particle_id))
    
    print("Finished in %s seconds" % (time.time() - start_time))
    
if __name__ == '__main__':
    correct_default_line_source_id()