#!/usr/bin/env python

"""
setup.py file for SWIG hoppet interface
"""

from distutils.core import setup, Extension


hoppet_v1_module = Extension('_hoppet_v1',
                           sources=['hoppet_v1_wrap.cxx'],
                           include_dirs = ['../src'],
                           swig_opts=['-c++'],
                           libraries=['hoppet_v1', 'gfortran'])

setup (name = 'hoppet_v1',
       version     = '1.3.0',
       author      = "Frederic Dreyer, Alexander Karlberg, Paolo Nason, Juan Rojo, Gavin P. Salam and Giulia Zanderighi",
       license     = 'GPLv3',
       description = """A python interface for Hoppet v1 created using SWIG""",
       ext_modules = [hoppet_v1_module],
       py_modules  = ["hoppet_v1"],
       )
