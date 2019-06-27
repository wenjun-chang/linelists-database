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
wavelengths = np.logspace(np.log10(0.3e-4), np.log10(30e-4), 4616)

plt.plot(wavelengths, database_absorption)
plt.plot(wavelengths, exocross_absorption)
plt.xscale('log')
plt.yscale('log')
plt.title('absorption cross section')
plt.grid(True)

plt.show()
