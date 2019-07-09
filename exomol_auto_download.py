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

############################

#fine
def get_molecule_links():
    url = r'http://exomol.com/data/molecules/'
    source_code = requests.get(url)
    plain_text = source_code.text
    soup = BeautifulSoup(plain_text, features='lxml')
    for link in soup.find_all('a', attrs={'class': 'list-group-item link-list-group-item'}):
        mol_name = link.get('href')
        href = 'http://exomol.com/data/molecules/' + mol_name
        print(mol_name, href)
        
        '''
        get_trans_files(href, mol_name)
        get_broad_files(href, mol_name)
        '''

#helper for getting trans files thx stackoverflow
def has_href_and_title_but_no_class(tag):
    return tag.has_attr('href') and tag.has_attr('title') and not tag.has_attr('class')
    

def get_trans_files(molecule_url, molecule_name):
    source_code = requests.get(molecule_url)
    plain_text = source_code.text
    soup = BeautifulSoup(plain_text, features='lxml')
    
    #find <a for all the links to the trans files for isotopologues of the molecules
    #at the main page for that molecule
    for link in soup.find_all('a', attrs={'class': 'list-group-item link-list-group-item'}):
        #print(link)
        href = molecule_url + '/' + link.get('href')
        #print(href)
        source_code = requests.get(href)
        plain_text = source_code.text
        soup = BeautifulSoup(plain_text, features='lxml')
        
        #go into the main page for data sets for the isotopologue
        for link in soup.find_all('a', attrs={'title': 'SAlTY'}):
            #print(link)
            href2 = href + '/' + link.get('href')
            #print(href2)
            source_code = requests.get(href2)
            plain_text = source_code.text
            soup = BeautifulSoup(plain_text, features='lxml')
            
            file_num = 1
            #ideally at the main page for the desired data type
            for link in soup.find_all(has_href_and_title_but_no_class):
                #1. trans
                #2. partition
                #print(link)
                if link is not None: 
                    href = link.get('href')
                    #the trans file
                    
                    if href.endswith('.trans.bz2'):
                        print(href)
                        trans_link = r'http://exomol.com' + href 
                        download_bz2_file(trans_link, molecule_name + '_trans_' + str(file_num))
                        file_num += 1
                    '''    
                    #the states file
                    elif href.endswith('.states.bz2'):
                        #print(href)
                        states_link = r'http://exomol.com' + href
                        download_bz2_file(states_link, molecule_name + '_states')
                        
                    #the partition file
                    elif href.endswith('.pf'):
                        #print(href)
                        partitions_link = r'http://exomol.com' + href
                        download_file(partitions_link, molecule_name + '_partitions')
                    '''
                    
#fine        
def get_broad_files(molecule_url, molecule_name):           
    source_code = requests.get(molecule_url)
    plain_text = source_code.text
    soup = BeautifulSoup(plain_text, features='lxml')   
   
    #get the links for the broadening parameter files and download them
    for link in soup.find_all(has_href_and_title_but_no_class):
        if link is not None:
            href = link.get('href')
            #H2 file
            if href.endswith('H2.broad'):
                H2_link = r'http://exomol.com' + href
                download_file(H2_link, molecule_name + '_H2_broad')
            #He file
            elif href.endswith('He.broad'):
                He_link = r'http://exomol.com' + href
                download_file(He_link, molecule_name + '_He_broad')
              
#fine
def download_file(url, outfile_name):
    #get file info
    file = request.urlopen(url)
    file_info = str(file.read())[2:-1]
    print(file_info)
    file_lines = file_info.split('\\n')[:-1]
    print(file_lines)
    #store file info in outfile
    outfile = open('{}.txt'.format(outfile_name), 'w')
    for line in file_lines:
        outfile.write(line + '\n')
    outfile.close()
    file.close()

#fine
def download_bz2_file(bz2_url, outfile_name):
    #get file info
    file = request.urlopen(bz2_url)
    CHUNK = 16 * 1024
    
    decompressor = bz2.BZ2Decompressor()
    with open('{}.txt'.format(outfile_name), 'wb') as outfile:
        while True:
            chunk = file.read(CHUNK)
            data = decompressor.decompress(chunk)
            print(data)
            if not chunk:
                break
            outfile.write(data)
    file.close()
           
    
##################       
   
def main():
    start_time = time.time()
    #get_molecule_links()
    get_trans_files(r'http://exomol.com/data/molecules/PH3', 'PH3')
    #get_broad_files(r'http://exomol.com/data/molecules/PH3', 'PH3')
    print("Finished in %s seconds" % (time.time() - start_time))
    
if __name__ == '__main__':
    main()