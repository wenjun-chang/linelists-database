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
import partition_calculator
import exomol_import

############################
#populate these after populating PH3
astro_molecules = ['H2O', 'CO', 'CH4', 'HCN', 'NH3', 'TiO', 'VO'] #'CO2'
#CO2 weird ass format lets only use zack and it is weird formatting need specify
#spedcial case CaO got 2 states files need hardcode

#fine
def get_molecule_links():
    url = r'http://exomol.com/data/molecules/'
    source_code = requests.get(url)
    plain_text = source_code.text
    soup = BeautifulSoup(plain_text, features='lxml')
    for link in soup.find_all('a', attrs={'class': 'list-group-item link-list-group-item'}):
        mol_name = link.get('href')
        href = 'http://exomol.com/data/molecules/' + mol_name
        #print(mol_name, href)
        
        if mol_name in astro_molecules:
            print(mol_name, href)
            #if no broad files, need to do something in exomol_import since broad files are not loadable
            #such as skipping loading broad file section and just assign all stuff to default
            has_broad_files = get_broad_files(href, mol_name)
            
            broad_H2_filepath = None
            broad_He_filepath = None
            if has_broad_files is True: 
                broad_H2_filepath = '/home/toma/Desktop/linelists-database/{}_H2_broad.txt{}'.format(mol_name)
                broad_He_filepath = '/home/toma/Desktop/linelists-database/{}_He_broad.txt{}'.format(mol_name)
            '''
            ######################
            #specially designed for CO2 weird format
            if mol_name == 'CO2':
                source_code = requests.get(r'http://exomol.com/data/molecules/CO2/')
                plain_text = source_code.text
                soup = BeautifulSoup(plain_text, features='lxml')
                
                for link in soup.find_all('a', attrs={'class': 'list-group-item link-list-group-item'}):
                    #print(link)
                    versions= []
                    iso_name = ''
                    atoms = link.get('href').split('-')
                    
                    for atom in atoms:
                        iso_name = '(' + atom + ')'  
                        href = molecule_url + '/' + link.get('href')
                        #print(href)
                        source_code = requests.get(href)
                        plain_text = source_code.text
                        soup = BeautifulSoup(plain_text, features='lxml')
                    
                        for link in soup.find_all('a', attrs={'title': 'Zak'}):
                            href2 = href + '/' + link.get('href')
                            #print(href2)
                            source_code = requests.get(href2)
                            plain_text = source_code.text
                            soup = BeautifulSoup(plain_text, features='lxml')
                        #i dont understand....
            ######################    
            '''
            suffixes = get_trans_files(href, mol_name)
            
            #use exomol_import file to import data into the database
            for i in range(len(suffixes)): #suffixes[i] = [iso_name, (version_name, trans_file_num, DEFAULT_GAMMA, DEFAULT_N)]
                iso_name = suffixes[i][0]
                versions_and_file_nums_and_default_and_ref_links = suffixes[i][1]
                version_name, trans_file_num, DEFAULT_GAMMA, DEFAULT_N, ref_links = versions_and_file_nums_and_default_and_ref_links
                
                states_filepath = '/home/toma/Desktop/linelists-database/' + mol_name + '_states_' + iso_name + '_' + version_name
                partitions_filepath = '/home/toma/Desktop/linelists-database/' + mol_name + '_partitions_' + iso_name + '_' + version_name
                trans_filepath_without_file_number = '/home/toma/Desktop/linelists-database/' + mol_name + '_trans_' + iso_name + '_' + version_name + '_'
                
                exomol_import.import_exomol_data(mol_name, iso_name, version_name, trans_filepath_without_file_number, \
                                                 states_filepath, partitions_filepath, broad_H2_filepath, broad_He_filepath, \
                                                 DEFAULT_GAMMA, DEFAULT_N, trans_file_num, ref_links)
            
        
#helper for getting trans files thx stackoverflow
def has_href_and_title_but_no_class(tag):
    return tag.has_attr('href') and tag.has_attr('title') and not tag.has_attr('class') 

def get_trans_files(molecule_url, molecule_name):
    suffixes = []
    
    source_code = requests.get(molecule_url)
    plain_text = source_code.text
    soup = BeautifulSoup(plain_text, features='lxml')
    
    #find <a for all the links to the main page for isotopologues of the molecules
    ###at the main page for that molecule listing all isotopologues
    for link in soup.find_all('a', attrs={'class': 'list-group-item link-list-group-item'}):
        #print(link)
        versions_and_file_nums_and_default_and_ref_links = []
        iso_name = ''
        atoms = link.get('href').split('-')
        for atom in atoms:
            iso_name = '(' + atom + ')'  
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
                href2 = href + '/' + link.get('href')
                #print(href2)
                source_code = requests.get(href2)
                plain_text = source_code.text
                soup = BeautifulSoup(plain_text, features='lxml')
                
                #check if files are properly present for each source
                for link in soup.find_all(has_href_and_title_but_no_class):
                    if link is not None: 
                        href = link.get('href')
                        
                        if href.endswith('.trans.bz2'):
                            has_trans = True
                        elif href.endswith('.states.bz2'):      
                            has_states = True
                        elif href.endswith('.pf'):      
                            has_partitions = True
                    
                if False in [has_states, has_trans]:
                    #skip this loop--don't download since there are insufficient file types
                    continue                
                
                file_num = 0
                DEFAULT_GAMMA = None
                DEFAULT_N = None
                #ideally at the main page for each source for each isotopologues
                #download all the trans, states, and partitions file
                for link in soup.find_all(has_href_and_title_but_no_class):
                    #1. trans
                    #2. partition
                    #3. states
                    #4. default gamma and N
                    #print(link)
                    if link is not None: 
                        href = link.get('href')
                        
                        #need to get def file and get default gamma and N value 
                        if href.endswith('.def'):
                            def_link = r'http://exomol.com' + href
                            print(def_link)
                            def_file_name = molecule_name + '_def_' + iso_name + '_' + version_name
                            download_file(def_link, def_file_name)
                            with open('/home/toma/Desktop/linelists-database/' + def_file_name + '.txt') as def_file: 
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
                            print(DEFAULT_GAMMA, DEFAULT_N)
                            
                        #the trans file
                        if href.endswith('.trans.bz2'):
                            file_num += 1
                            trans_link = r'http://exomol.com' + href
                            print(trans_link)
                            download_bz2_file(trans_link, molecule_name + '_trans_' + iso_name + '_' + version_name + '_' + str(file_num))
                        '''
                        #the states file
                        elif href.endswith('.states.bz2'):
                            states_link = r'http://exomol.com' + href
                            print(states_link)
                            download_bz2_file(states_link, molecule_name + '_states_' + iso_name + '_' + version_name)
                            
                        #the partition file
                        elif href.endswith('.pf'):
                            partitions_link = r'http://exomol.com' + href
                            print(partitions_link)
                            download_file(partitions_link, molecule_name + '_partitions_' + iso_name + '_' + version_name)
                        '''
                        
                #if no partition file...can only use this instead of downloading the partition file form exomol in the future
                if has_partitions is False: 
                    #calculate partiton function and write to a file
                    states_filepath = '/home/toma/Desktop/linelists-database/' + molecule_name + '_states_' + iso_name + '_' + version_name
                    partition_filepath = '/home/toma/Desktop/linelists-database/' + molecule_name + '_partitions_' + iso_name + '_' + version_name
                    #this fucntion will calculate partition and save it to the same filepath as if the partition file exists
                    partition_calculator.calculate_partition(states_filepath, partition_filepath, 10000)
                
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
                            seen.append(link)
                for link in seen: 
                    ref_links += '%s||' % link
                ref_links = ref_links[:-2]
                
                #append the fucntioning version and number of trans files to returning output 
                versions_and_file_nums_and_default_and_ref_links.append((version_name, file_num, DEFAULT_GAMMA, DEFAULT_N, ref_links))
            
        '''
        #if number of partition files does not match number of useful version available
        if num_partitions != len(versions):
            #do sth
            pass
        if num_partitions == 0:
            #delete all files for this isotopologue if no single partition file exist for this isotopologues
            for filename in os.listdir('/home/toma/Desktop/linelists-database'):
                if iso_name in filename:
                    os.remove('/home/toma/Desktop/linelists-database/' + filename)
        '''
        #if not at least one of parition, states, and trans file is existing then delete all
        #if only one partition file, chnage partition filename to general name
        
        versions = [i[0] for i in versions_and_file_nums_and_default_and_ref_links]
        print('version for', iso_name, 'includes', *versions)
        suffixes.append([iso_name, versions_and_file_nums_and_default_and_ref_links])
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
    file = request.urlopen(url)
    file_info = file.read()
    #store file info in outfile
    outfile = open('{}.txt'.format(outfile_name), 'wb')
    outfile.write(file_info)
    outfile.close()
    file.close()

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
    os.system('bzip2 -d /home/toma/Desktop/linelists-database/{}'.format(compressed))
    
    '''
    #can run out of memory
    file = request.urlopen(bz2_url)
    #decompressor = bz2.BZ2Decompressor()
    compressed_data = file.read()
    data = bz2.decompress(compressed_data)
    #file_lines = data.decode()
    #store file info in outfile
    outfile = open('{}.txt'.format(outfile_name), 'wb')
    outfile.write(data)
    outfile.close()
    '''
    
##################       
   
def main():
    start_time = time.time()
    #get_molecule_links()
    #get_trans_files(r'http://exomol.com/data/molecules/PH3', 'PH3')
    #get_broad_files(r'http://exomol.com/data/molecules/PH3', 'PH3')
    #download_bz2_file(r'http://exomol.com/db/PH3/31P-1H3/SAlTY/31P-1H3__SAlTY__00000-00100.trans.bz2', 'test')
    print("Finished in %s seconds" % (time.time() - start_time))
    
if __name__ == '__main__':
    main()