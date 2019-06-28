#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 16:40:46 2019

@author: toma
"""

import numpy as np
import matplotlib.pyplot as plt

exocross_absorption = np.load('/home/toma/Desktop/exocross_sigmas.npy')
database_absorption = np.load('/home/toma/Desktop/absorption.npy')
wavelengths = np.logspace(np.log10(300), np.log10(3000), 4616)

plt.plot(wavelengths, database_absorption)
plt.plot(wavelengths, exocross_absorption)
plt.xscale('log')
plt.yscale('log')
plt.title("computed at 2000K, 0.1 atm for '(12C)(16O)' using HITRAN 2016 lines at the range 300 um to 3000 um using logspace display", fontsize='x-large')
plt.xlabel('wavelength (um)', fontsize='x-large')
plt.ylabel('absorption cross section (cm^2)', fontsize='x-large')
plt.get_xaxis().get_major_formatter().set_useOffset(False)
plt.grid(True)

plt.show()
