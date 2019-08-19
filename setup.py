#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 11:43:54 2019

@author: toma
"""

#setup file for cython optimization 

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
#import Cython.Compiler.Options
#directive_defaults = Cython.Compiler.Options.get_directive_defaults()

#directive_defaults['linetrace'] = True
#directive_defaults['binding'] = True

'''
setup(
    ext_modules = cythonize("compute_absorption_cython.pyx")
)
'''

ext_modules=[ Extension("compute_absorption_c",
              ["compute_absorption_cython.pyx", "Faddeeva.cc"],
              libraries=["m"],
              extra_compile_args = ["-ffast-math"], 
              language="c++", 
              define_macros=[('CYTHON_TRACE', '1')])]

setup(
      name = "compute_absorption_c",
      cmdclass = {"build_ext": build_ext},
      compiler_directives={'linetrace':True, 'binding':True},
      ext_modules = ext_modules
      )

