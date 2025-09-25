#!/usr/bin/env python

"""
Legacy setup.py file for hoppet interface.

Note: This file is kept for backwards compatibility but is no longer used.
The build system now uses pyproject.toml with scikit-build-core and pybind11.
"""

# This file is legacy and not used in the current build system.
# The actual Python interface is built using:
# - pyproject.toml (build configuration)
# - CMakeLists.txt (build system)
# - hoppet_pybind11.cpp (pybind11 bindings)

print("WARNING: This setup.py file is legacy and not used.")
print("The Python interface is now built using pyproject.toml with scikit-build-core and pybind11.")
print("Please use 'pip install .' from the root directory to build the package.")
