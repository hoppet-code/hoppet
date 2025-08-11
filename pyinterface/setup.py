#!/usr/bin/env python

"""
setup.py file for SWIG hoppet interface
"""

from distutils.core import setup, Extension


hoppet_module = Extension('_hoppet',
                           sources=['hoppet_wrap.cxx'],
                           include_dirs = ['../src','.'],
                           swig_opts=['-c++'],
                           libraries=['hoppet', 'gfortran'])

setup (name = 'hoppet',
       version     = '2.0.0',
       author      = "Frederic Dreyer, Alexander Karlberg, Paolo Nason, Juan Rojo, Gavin P. Salam and Giulia Zanderighi",
       author_email= "gavin.salam@physics.ox.ac.uk",
       url         = "https://github.com/hoppet-code/hoppet",
       license     = 'GPLv3',
       description = """A python interface for Hoppet v1 created using SWIG""",
       ext_modules = [hoppet_module],
       py_modules  = ["hoppet"],
       )
