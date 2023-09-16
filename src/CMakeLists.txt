## Experimental CMakeLists.txt file for Hoppet (library only, installation still rocky)
## 
## mkdir build; cd build
## cmake [-DCMAKE_INSTALL_PREFIX=/opt/myproject] [-DUSER_EXACT_COEFF=ON] ..
## make -j  [or cmake --build -j .
## cmake --install .
##

cmake_minimum_required(VERSION 3.10)
# set the project name
project(Hoppet VERSION 1.2.1)
# specify that we're using fortran
enable_language(Fortran)
set(CMAKE_Fortran_FLAGS "-O3")

# add an option to decide whether to use the exact coefficient functions
option(USE_EXACT_COEF "Use exact coefficient functions" OFF)

# create a list of all .f and .f90 files in the current directory
file(GLOB SOURCES *.f90 *.f)

# alternative is a manual list (NB: this specific list is not complete)
# set(SOURCES
#   hplog.f xpij2e.f xpns2e.f assertions.f90 coefficient_functions.f90
#   convolution.f90 dglap_choices.f90 dglap_holders.f90 dglap_objects.f90
#   dummy_pdfs.f90 evolution.f90 hoppet_v1.f90 integrator.f90
#   interpolation.f90 new_as.f90 pdf_general.f90 pdf_representation.f90
#   pdf_tabulate.f90 qcd.f90 qcd_coupling.f90 random.f90 runge_kutta.f90
#   sort.f90 special_functions.f90 splitting_functions.f90
#   splitting_functions_nnlo.f90 streamlined_interface.f90 types.f90
#   warnings_and_errors.f90 welcome_message.f90 xa2hgp.f90 xpij2n.f90
#   xpij2p.f90 xpns2n.f90 xpns2p.f90)

# and additionally in the directories used for conditional compilation
file(GLOB EXACT_COEF_SOURCES exact-coefs/*.f exact-coefs/*.f90)
file(GLOB NO_EXACT_COEF_SOURCES no-exact-coefs/*.f no-exact-coefs/*.f90)

if(<variable-or-expression>)
  # commands to execute if the condition is true
endif()

# if USE_EXACT_COEF is ON, use the exact coefficient functions sources
if(USE_EXACT_COEF)
  set(SOURCES ${SOURCES} ${EXACT_COEF_SOURCES})
else()
  set(SOURCES ${SOURCES} ${NO_EXACT_COEF_SOURCES})
endif()

# 
# # build a library from those sources
add_library(hoppet SHARED ${SOURCES})


# arrange for installation of the library
# NB: installation of modules is not yet clear
install(TARGETS hoppet DESTINATION lib)
install(FILES hoppet_v1.h DESTINATION include)


#set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod_files)
#install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include)

#install(TARGETS hoppet PUBLIC_HEADER DESTINATION lib)
# #-------------------------------------------------------------------------
# # main code 
# add_library(hoppet SHARED
#   hplog.f xpij2e.f xpns2e.f assertions.f90 coefficient_functions.f90
#   convolution.f90 dglap_choices.f90 dglap_holders.f90 dglap_objects.f90
#   dummy_pdfs.f90 evolution.f90 hoppet_v1.f90 integrator.f90
#   interpolation.f90 new_as.f90 pdf_general.f90 pdf_representation.f90
#   pdf_tabulate.f90 qcd.f90 qcd_coupling.f90 random.f90 runge_kutta.f90
#   sort.f90 special_functions.f90 splitting_functions.f90
#   splitting_functions_nnlo.f90 streamlined_interface.f90 types.f90
#   warnings_and_errors.f90 welcome_message.f90 xa2hgp.f90 xpij2n.f90
#   xpij2p.f90 xpns2n.f90 xpns2p.f90)
