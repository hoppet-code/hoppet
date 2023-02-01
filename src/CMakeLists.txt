cmake_minimum_required(VERSION 3.10)
# set the project name
project(Hoppet VERSION 1.2.1)
# specify that we're using fortran
enable_language(Fortran)
set(CMAKE_Fortran_FLAGS "-O3")
#-------------------------------------------------------------------------
# main code 
add_library(hoppet SHARED
  hplog.f xpij2e.f xpns2e.f assertions.f90 coefficient_functions.f90
  convolution.f90 dglap_choices.f90 dglap_holders.f90 dglap_objects.f90
  dummy_pdfs.f90 evolution.f90 hoppet_v1.f90 integrator.f90
  interpolation.f90 new_as.f90 pdf_general.f90 pdf_representation.f90
  pdf_tabulate.f90 qcd.f90 qcd_coupling.f90 random.f90 runge_kutta.f90
  sort.f90 special_functions.f90 splitting_functions.f90
  splitting_functions_nnlo.f90 streamlined_interface.f90 types.f90
  warnings_and_errors.f90 welcome_message.f90 xa2hgp.f90 xpij2n.f90
  xpij2p.f90 xpns2n.f90 xpns2p.f90)
