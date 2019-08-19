#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 16:40:46 2019

@author: toma
"""

import numpy as np
import matplotlib.pyplot as plt

#this file contains the fucntion that plots absorption cross section against wavelength
#it also can compare exocross and database absorption on the same plot

def plot_absorption(wavenum_file, absorption_file): 
    wavenums = np.loadtxt(wavenum_file, usecols=0, unpack=True)
    wavelengths = (1 / wavenums) * 1e4
    absorptions = np.loadtxt(absorption_file, usecols=0, unpack=True)
    
    plt.xscale('log')
    plt.yscale('log')
    plt.plot(wavelengths, absorptions, label='Absorption Plot')
    
    plt.gca().legend(('database absorption',), fontsize='xx-large')
    plt.title("Absorption cross section plot", fontsize='xx-large')
    plt.xlabel('Wavelength (um)', fontsize='xx-large')
    plt.ylabel('Absorption Cross Section (cm^2)', fontsize='xx-large')
    plt.grid(True)
    plt.show()
    
if __name__ == '__main__': 
    #database_absorption = np.load('/home/toma/Desktop/absorption.npy')
    database_absorption = np.load('/home/toma/Desktop/telluric_spectrum_model.npy')
    wavelengths, telluric_spectrum = np.loadtxt('/home/toma/Desktop/Telluric Spectrum', skiprows=1, usecols=(2, 3), unpack=True)
    #wavelengths = (1 / wavenums) * 1e4
    
    #plt.subplot(1, 2, 1)
    database_model = plt.plot(wavelengths, database_absorption, label='database model')
    #plt.subplot(1, 2, 2)
    telluric_spectrum = plt.plot(wavelengths, telluric_spectrum, label='telluric spectrum')
    
    #plt.plot(wavelengths, database_absorption/telluric_spectrum)
    plt.xscale('log')
    plt.yscale('log')
    plt.gca().legend(('database_model', 'telluric_spectrum'), fontsize='xx-large')
    #plt.title("Telluric Spectrum modeled at 270K, 0.5 atm using HITRAN 2016 (1H)2(16O) line list of <10% error", fontsize='xx-large')
    plt.xlabel('Wavelength (um)', fontsize='xx-large')
    plt.ylabel('Spectrum (e)', fontsize='xx-large')
    plt.grid(True)
    
    plt.show()