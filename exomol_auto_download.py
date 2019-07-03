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

#need write more, even helper function for open url to save lines
def get_trans_files(molecule_url, molecule_name):
    source_code = requests.get(molecule_url)
    plain_text = source_code.text
    soup = BeautifulSoup(plain_text, features='lxml')
    
    #find <a for all the links to the trans files for isotopologues of the molecules
    for link in soup.find_all('a', attrs={'class': 'list-group-item link-list-group-item'}):
        print(link)
        href = molecule_url + link.get('href')
        source_code = requests.get(href)
        plain_text = source_code.text
        soup = BeautifulSoup(plain_text, features='lxml')
        for link in soup.find_all('a', attrs={'title': 'SAlTY'}):
            ####more links afterwards
            file_link = molecule_url + link.get('href')
            download_file(file_link, molecule_name + '_transitions')
            pass
#fine, file a bit messy though          
def get_broad_files(molecule_url, molecule_name):            
    source_code = requests.get(molecule_url)
    plain_text = source_code.text
    soup = BeautifulSoup(plain_text, features='lxml')    

    #get the links for the broadening parameter files and download them
    for link in soup.find_all('a', attrs={'href': '/db/PH3/31P-1H3/31P-1H3__H2.broad'}):
        H2_link = r'http://exomol.com' + link.get('href')        
        print(H2_link)
        download_file(H2_link, molecule_name + '_H2_broad')
        
    for link in soup.find_all('a', attrs={'href': '/db/PH3/31P-1H3/31P-1H3__He.broad'}):
        He_link = r'http://exomol.com' + link.get('href')
        print(He_link)
        download_file(He_link, molecule_name + '_He_broad')
        
#fine
def download_file(url, outfile_name):
    #get file info
    file = request.urlopen(url)
    file_info = str(file.read())
    file_lines = file_info.split('\\n')
    
    #store file info in outfile
    outfile = open('{}.txt'.format(outfile_name), 'w')
    for line in file_lines:
        outfile.write(line + '\n')
    outfile.close()
    
    
##################        
def main():
    
    start_time = time.time()
    
    #get_molecule_links()
    get_broad_files(r'http://exomol.com/data/molecules/PH3', 'PH3')
    
    print("Finished in %s seconds" % (time.time() - start_time))
    
if __name__ == '__main__':
    main()


