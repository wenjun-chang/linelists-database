#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 16:40:46 2019

@author: toma
"""

import numpy as np
import matplotlib.pyplot as plt

exocross_absorption = np.load('/home/toma/Desktop/sigmas_CO_1000K_0.1bar.npy')
database_absorption = np.load('/home/toma/Desktop/absorption.npy')
#wavelengths = np.logspace(np.log10(300), np.log10(3000), 4616)
wavelengths = np.load('/home/toma/Desktop/wavelengths.npy')

#plt.subplot(1, 2, 1)
plt.plot(wavelengths, database_absorption)
#plt.subplot(1, 2, 2)
plt.plot(wavelengths, exocross_absorption)
plt.xscale('log')
plt.yscale('log')
plt.title("computed at 1000K, 0.1 atm for '(12C)(16O)' using ExoMol Li2015 lines at the range 0.943 um to 2.44 um using logspace display", fontsize='medium')
plt.xlabel('wavelength (um)', fontsize='medium')
plt.ylabel('absorption cross section (cm^2)', fontsize='medium')
plt.grid(True)

plt.show()
