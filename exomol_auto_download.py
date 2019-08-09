#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue Jul  2 15:30:48 2019
@author: toma
"""

import numpy as np
import time
from urllib import request
import requests
from bs4 import BeautifulSoup
import bz2
import os
from partition_calculator import calculate_partition
from exomol_import import import_exomol_data
from query_functions import sql_order

############################
astro_molecules = ['HCN', 'CH4', 'CO', 'VO', 'TiO', 'NH3', 'H2O', 'CO2'] 
#'CO2' has some problem,'HCN', 'CH4', 'CO', 'VO', 'TiO', 'NH3', 'H2O' (without superline) done
#spedcial case CaO got 2 states files need hardcode

#fine
def populate_all_exomol():
    start_time = time.time()
    
    url = r'http://exomol.com/data/molecules/'
    source_code = requests.get(url)
    plain_text = source_code.text
    soup = BeautifulSoup(plain_text, features='lxml')
    for link in soup.find_all('a', attrs={'class': 'list-group-item link-list-group-item'}):
        mol_name = link.get('href')
        href = 'http://exomol.com/data/molecules/' + mol_name
        #print(mol_name, href)
        
        if mol_name not in astro_molecules and mol_name != 'PH3' and mol_name != 'C2H4':
            print(mol_name, href)
            #if no broad files, need to do something in exomol_import since broad files are not loadable
            #such as skipping loading broad file section and just assign all stuff to default
            has_broad_files = get_broad_files(href, mol_name)
            
            broad_H2_filepath = None
            broad_He_filepath = None
            if has_broad_files is True: 
                broad_H2_filepath = '/home/toma/Desktop/linelists-database/{}_H2_broad'.format(mol_name)
                broad_He_filepath = '/home/toma/Desktop/linelists-database/{}_He_broad'.format(mol_name)
            
            suffixes = get_trans_files(href, mol_name)
            
            #use exomol_import file to import data into the database
            for i in range(len(suffixes)): #suffixes[i] = [iso_name, [(version_name, trans_file_num, DEFAULT_GAMMA, DEFAULT_N), ..., ...]]
                iso_name = suffixes[i][0]
                if suffixes[i][1] == []:
                    print('sdsadasbfabfadbvbdvasvbisavbisfv')
                    continue
                versions_and_file_nums_and_default_and_ref_links_list = suffixes[i][1]
                num_versions = len(versions_and_file_nums_and_default_and_ref_links_list)
                for i in range(num_versions):  
                    version_name, trans_file_num, DEFAULT_GAMMA, DEFAULT_N, ref_links = versions_and_file_nums_and_default_and_ref_links_list[i]
                    
                    states_filepath = '/home/toma/Desktop/linelists-database/' + mol_name + '_states_' + iso_name + '_' + version_name
                    partitions_filepath = '/home/toma/Desktop/linelists-database/' + mol_name + '_partitions_' + iso_name + '_' + version_name
                    trans_filepath_without_file_number = '/home/toma/Desktop/linelists-database/' + mol_name + '_trans_' + iso_name + '_' + version_name + '_'
                    
                    import_exomol_data(mol_name, iso_name, version_name, trans_filepath_without_file_number, states_filepath, \
                                       partitions_filepath, broad_H2_filepath, broad_He_filepath, DEFAULT_GAMMA, DEFAULT_N, \
                                       trans_file_num, ref_links)
               
    print("Finished in %s seconds" % (time.time() - start_time))
    
#helper for getting trans files thx stackoverflow
def has_href_and_title_but_no_class(tag):
    return tag.has_attr('href') and tag.has_attr('title') and not tag.has_attr('class')

#helper for getting trans files thx stackoverflow
def has_href_and_class_but_no_title(tag):
    return tag.has_attr('href') and tag.has_attr('class') and not tag.has_attr('title') 

def get_trans_files(molecule_url, molecule_name):
    suffixes = []
    
    source_code = requests.get(molecule_url)
    plain_text = source_code.text
    soup = BeautifulSoup(plain_text, features='lxml')
    
    DEFAULT_GAMMA = None
    DEFAULT_N = None
    #find <a for all the links to the main page for isotopologues of the molecules
    ###at the main page for that molecule listing all isotopologues
    for link in soup.find_all('a', attrs={'class': 'list-group-item link-list-group-item'}):
        #print(link)
        versions_and_file_nums_and_default_and_ref_links = []
        iso_name = ''
        atoms = link.get('href').split('-')
        for atom in atoms:
            if atom[-1:].isdigit() is True: #if atom is like 16O2
                num_atom = atom[-1:]
                atom = atom[:-1]
                iso_name += '(' + atom + ')' + num_atom
            else: #if atom is like 12C
                iso_name += '(' + atom + ')'
        print(iso_name)
        href = molecule_url + '/' + link.get('href')
        #print(href)
        source_code = requests.get(href)
        plain_text = source_code.text
        soup = BeautifulSoup(plain_text, features='lxml')
        
        #find <a for all the links to the main page for different versions of the isotopologue
        ###at main page for each isotopologue
        ###############slight modification needed for full automation
        for link in soup.find_all('a', attrs={'class': 'list-group-item link-list-group-item'}):
            version_name = link.get('title')
                        
            has_states = False
            has_trans = False
            has_partitions = False
                            
            #make sure it is not an external link
            if version_name.startswith('xsec') is False:
                #print(link)
                version_name = 'EXOMOL_' + version_name
                href2 = ''
                href2 = href + '/' + link.get('href')
                print(href2)
                source_code = requests.get(href2)
                plain_text = source_code.text
                soup = BeautifulSoup(plain_text, features='lxml')
                
                #check if files are properly present for each source
                for link in soup.find_all(has_href_and_title_but_no_class):
                    if link is not None: 
                        href3 = link.get('href')
                        
                        
                        if href3.endswith('.trans.bz2'):
                            has_trans = True
                        elif href3.endswith('.states.bz2'):      
                            has_states = True
                        elif href3.endswith('.pf'):      
                            has_partitions = True
                    
                if False in [has_states, has_trans]:
                    #skip this loop--don't download since there are insufficient file types
                    print('wooooooooooooooooooooooooooooooo')
                    continue                
                
                file_num = 0
                #ideally at the main page for each source for each isotopologues
                #download all the trans, states, and partitions file
                #need to get def file and get default gamma and N value 
                
                #get default values from def file
                for link in soup.find_all(has_href_and_class_but_no_title):
                    if link is not None: 
                        href4 = link.get('href')
                        if href4.endswith('.def'):
                            def_link = r'http://exomol.com' + href4
                            #print(def_link)
                            def_file_name = molecule_name + '_def_' + iso_name + '_' + version_name
                            has_def = download_file(def_link, def_file_name)
                            if has_def is True: 
                                with open('/home/toma/Desktop/linelists-database/' + def_file_name) as def_file: 
                                    for line in def_file: 
                                        data = line.strip().split()
                                        if 'Default' in data: 
                                            if 'Lorentzian' in data:                 
                                                DEFAULT_GAMMA = float(data[0])
                                            else:
                                                DEFAULT_N = float(data[0])
                                def_file.close()
                                #delete def file cuz no longer needed
                                os.remove('/home/toma/Desktop/linelists-database/' + def_file_name)
                                #print(DEFAULT_GAMMA, DEFAULT_N)
                            else: #when the stupid exomol url is invalid
                                DEFAULT_GAMMA = DEFAULT_GAMMA
                                DEFAULT_N = DEFAULT_N
                
                for link in soup.find_all(has_href_and_title_but_no_class):
                    #1. trans
                    #2. partition
                    #3. states
                    #print(link)
                    downloaded = ['HCN', 'CH4', 'CO', 'VO', 'TiO', 'NH3', 'H2O', 'CO2'] #########
                    if link is not None: 
                        href5 = link.get('href')                            
                        #the trans file
                        if href5.endswith('.trans.bz2'):
                            file_num += 1
                            if molecule_name not in downloaded: #####
                                trans_link = r'http://exomol.com' + href5
                                #print(trans_link)
                                download_bz2_file(trans_link, molecule_name + '_trans_' + iso_name + '_' + version_name + '_' + str(file_num))
                        
                        #the states file
                        elif href5.endswith('.states.bz2'):
                            if molecule_name not in downloaded: #####
                                states_link = r'http://exomol.com' + href5
                                #print(states_link)
                                if molecule_name == 'CaO': #specially hardcoded for CaO which has 2 states files
                                    if href5.endswith('.abinitio.states.bz2'):
                                        continue #dont download the theoratical states file
                                download_bz2_file(states_link, molecule_name + '_states_' + iso_name + '_' + version_name)
                            
                        #the partition file
                        elif href5.endswith('.pf'):
                            if molecule_name not in downloaded: #####
                                has_partitions = True
                                partitions_link = r'http://exomol.com' + href5
                                #print(partitions_link)
                                download_file(partitions_link, molecule_name + '_partitions_' + iso_name + '_' + version_name)                
                                    
                #if no partition file...can only use this instead of downloading the partition file form exomol in the future
                if has_partitions is False: 
                    if molecule_name not in downloaded: #####
                        #calculate partiton function and write to a file
                        states_filepath = '/home/toma/Desktop/linelists-database/' + molecule_name + '_states_' + iso_name + '_' + version_name
                        partition_filepath = '/home/toma/Desktop/linelists-database/' + molecule_name + '_partitions_' + iso_name + '_' + version_name
                        #this fucntion will calculate partition and save it to the same filepath as if the partition file exists
                        calculate_partition(states_filepath, partition_filepath, 10000)
                
                #find the link to the reference paper...now just fetch all the links
                #still dont know how to handle this
                #:::<span class="noprint"> [<a href="http://dx.doi.org/10.1093/mnras/stu2246">link to article</a>]</span>
                #reference links will be separated by ||
                #i.e. link1 || link2 ||...||linkn-1||linkn
                ref_links = ''
                seen = []
                for link in soup.find_all('span', attrs={'class': 'noprint'}):
                    if link is not None:
                        link = str(link)
                        link = link[link.find('[') + 1 : link.find(']')]
                        link = link[link.find('"') + 1 :]
                        link = link[: link.find('"')]
                        if link not in seen: 
                            #print(link)
                            seen.append(link)
                for link in seen: 
                    ref_links += '%s||' % link
                ref_links = ref_links[:-2]
                
                #append the fucntioning version and number of trans files to returning output 
                versions_and_file_nums_and_default_and_ref_links.append((version_name, file_num, DEFAULT_GAMMA, DEFAULT_N, ref_links))
            
        print(iso_name, versions_and_file_nums_and_default_and_ref_links)
        versions = [i[0] for i in versions_and_file_nums_and_default_and_ref_links]
        print('Versions for', molecule_name, 'includes', *versions)
        suffixes.append([iso_name, versions_and_file_nums_and_default_and_ref_links])
    print(suffixes)
    return suffixes #use this for full automation connect to exomol_import
        
#fine
#if no broad maybe go into get trans and find the header??
def get_broad_files(molecule_url, molecule_name):
    has_broad_files = False
         
    source_code = requests.get(molecule_url)
    plain_text = source_code.text
    soup = BeautifulSoup(plain_text, features='lxml')   
   
    #get the links for the broadening parameter files and download them
    for link in soup.find_all(has_href_and_title_but_no_class):
        if link is not None:
            href = link.get('href')
            #H2 file
            if href.endswith('H2.broad'):
                has_broad_files = True
                H2_link = r'http://exomol.com' + href
                download_file(H2_link, molecule_name + '_H2_broad')
            #He file
            elif href.endswith('He.broad'):
                has_broad_files = True
                He_link = r'http://exomol.com' + href
                download_file(He_link, molecule_name + '_He_broad')
    return has_broad_files
              
#fine
def download_file(url, outfile_name):
    #get file info
    try: 
        file = request.urlopen(url)
        file_info = file.read()
        #store file info in outfile
        outfile = open('{}'.format(outfile_name), 'wb')
        outfile.write(file_info)
        outfile.close()
        file.close()
        return True
    except Exception: 
        return False
    
#fine
def download_bz2_file(bz2_url, outfile_name):
    #get file info
    file = request.urlopen(bz2_url)
    file_info = file.read()
    #store file info in outfile
    compressed = '{}.bz2'.format(outfile_name)
    outfile = open(compressed, 'wb')
    outfile.write(file_info)
    outfile.close()
    file.close()
    os.system("bzip2 -d '/home/toma/Desktop/linelists-database/{}'".format(compressed))
    
    '''
    #can run out of memory
    file = request.urlopen(bz2_url)
    #decompressor = bz2.BZ2Decompressor()
    compressed_data = file.read()
    data = bz2.decompress(compressed_data)
    #file_lines = data.decode()
    #store file info in outfile
    outfile = open('{}'.format(outfile_name), 'wb')
    outfile.write(data)
    outfile.close()
    '''
    
##################   
if __name__ == '__main__':
    
    #disable autocommit to improve performance
    sql_order('SET autocommit = 0')
    sql_order('SET unique_checks = 0')
    sql_order('SET foreign_key_checks = 0')
    sql_order('SET sql_log_bin = 0')
    
    populate_all_exomol()
    
    #turn them back on
    sql_order('SET unique_checks = 1')
    sql_order('SET foreign_key_checks = 1')
    sql_order('SET sql_log_bin = 1')

    
