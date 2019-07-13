#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 16:40:46 2019

@author: toma
"""

import numpy as np
import matplotlib.pyplot as plt

temp = np.load('/home/toma/Desktop/temp.npy')
temp_exo = np.loadtxt('/home/toma/Desktop/only_one_line.xsec', usecols=1, unpack=True)

database_absorption = np.load('/home/toma/Desktop/absorption.npy')
wavenums, exocross_absorption = np.loadtxt('/home/toma/Desktop/output_1000_1d-1.xsec', unpack=True)
wavelengths = (1 / wavenums) * 1e4


#plt.subplot(1, 2, 1)
plt.plot(wavelengths, database_absorption/0.98654)
#plt.subplot(1, 2, 2)
plt.plot(wavelengths, exocross_absorption)


#plt.plot(wavelengths, temp/0.98654)
#plt.plot(wavelengths, temp_exo)

#ratio = exocross_absorption / (database_absorption/0.98654)
#plt.plot(wavelengths, np.abs(ratio))

plt.xscale('log')
plt.yscale('log')
plt.title("computed at 1000K, 0.1 atm for '(12C)(16O)' using ExoMol Li2015 lines at the range 0.943 um to 2.44 um using logspace display", fontsize='medium')
plt.xlabel('wavelength (um)', fontsize='medium')
plt.ylabel('absorption cross section (cm^2)', fontsize='medium')
plt.grid(True)

plt.show()
