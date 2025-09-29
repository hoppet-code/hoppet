# Usage instructions:
#
## mkdir build
## cmake -S . -B build [-DCMAKE_INSTALL_PREFIX=/opt/myproject] [-DHOPPET_USE_EXACT_COEF=ON] ..
## cmake --build  build -j 
## cmake --install build
## ctest --test-dir build [-j8]
#
# To see the HOPPET-specific CMake options, run `cmake -L ..`


cmake_minimum_required(VERSION 3.15)

# set the project name: note should be updated also in src/welcome_message.f90 and the documentation
project(hoppet VERSION 2.0.0 LANGUAGES C CXX Fortran)
set(PROJECT_VERSION_PRERELEASE "-devel")

set(HOPPET_VERSION "${PROJECT_VERSION}${PROJECT_VERSION_PRERELEASE}")
message(STATUS "HOPPET: HOPPET_VERSION=${HOPPET_VERSION}")

set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules ${CMAKE_MODULE_PATH})
include("GNUInstallDirs")

# set the C++ standard to c++14 for the example
set(CMAKE_CXX_STANDARD 14)

# add check to see if src/Makefile exists, and any src/*.mod, or any src/*.a
# If so, then bail out because the user probably tried a build with the old build system
file(GLOB MOD_FILES "${PROJECT_SOURCE_DIR}/src/*.mod")
file(GLOB A_FILES "${PROJECT_SOURCE_DIR}/src/*.a")
if(EXISTS "${PROJECT_SOURCE_DIR}/src/Makefile"
   OR (NOT "${MOD_FILES}" STREQUAL "")
   OR (NOT "${A_FILES}" STREQUAL ""))
  message(FATAL_ERROR "HOPPET: Old build system files detected (including potentially ${MOD_FILES} ${A_FILES})."
   " Please remove extraneous files or start from a fresh clone (you can also try 'make distclean' from the top-level directory, but we do not guarantee this will work)")
endif()

# Add options ==========================================================
option(HOPPET_USE_EXACT_COEF    "Use exact coefficient functions"  OFF)
option(HOPPET_BUILD_EXAMPLES    "Build examples" ON)
option(HOPPET_ENABLE_TESTING    "Enable testing. Requires building the examples." ON)
option(HOPPET_BUILD_BENCHMARK   "Build benchmark." ON)
option(HOPPET_BUILD_PYINTERFACE "Build Python interface." OFF)
option(HOPPET_BUILD_PYWHEELS    "Build Python wheels. (developer option)" OFF)
option(HOPPET_ENABLE_FPES       "Enable trapping for the floating point exceptions." OFF)
option(HOPPET_BUILD_TESTDIR     "Enable build of the code in the test directory." OFF)

# compile-time debugging info
option(HOPPET_ENABLE_DEBUG "Turn on debug compiler information [default=ON]" ON)


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RELEASE")
  #set(CMAKE_BUILD_TYPE "RelWithDebInfo")
message(STATUS "HOPPET: Build type not specified, use: ${CMAKE_BUILD_TYPE}")
endif()

if(HOPPET_BUILD_PYINTERFACE)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
elseif(HOPPET_BUILD_PYWHEELS)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
  #set(SKBUILD_PROJECT_VERSION "${HOPPET_VERSION}")
  #message(STATUS "Setting Python package version to: ${HOPPET_VERSION}")
endif()


string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE)
message(STATUS "HOPPET: CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")

set(HOPPET_PYTHON_PACKAGE_NAME "hoppet" CACHE STRING "The name of the installed hoppet python package.")
set(HOPPET_USE_PYTHON_SITEARCH OFF CACHE BOOL "Force the usage of ${Python_SITEARCH} if CMAKE_INSTALL_PREFIX was defined.")
set(HOPPET_CUSTOM_PYTHON_INSTALL "" CACHE STRING "Specify a custom python install path, can be absolute or relative to lib.")

if (HOPPET_ENABLE_DEBUG)
  add_compile_options("-g")
endif()


# Flags for compilers ==============================================
# GNU
#add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fstack-arrays>") # stack allocation of arrays seems to have no impact
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-ffree-line-length-none>") # allow lines>132 characters (as in recent std)
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-Wall>")
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-Wextra;-Wno-unused-parameter;-Wno-unused-variable>")
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-Wno-unused-dummy-argument;-Wno-maybe-uninitialized>")
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-Wno-unused;-Wno-uninitialized;-Wno-conversion>")
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-Wno-compare-reals;-Wno-function-elimination>")
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-Wno-integer-division>")
if (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS_DEBUG  "${CMAKE_Fortran_FLAGS_DEBUG} -g -fbacktrace -fcheck=all")
  # in qed-evolution, avoid a spurious target-lifetime warning by turning off that
  # warning just for that file
  set_source_files_properties(
    src/qed_evolution.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-target-lifetime"
  )
endif()

# Intel
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-diag-disable 8291;-diag-disable 10448;-diag-disable 8889;-diag-disable 6717;-diag-disable 8291;-diag-disable 7712;-diag-disable 5194>")
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>:-warn;-diag-disable 8291;-diag-disable 8889;-diag-disable 6717;-diag-disable 8291;-diag-disable 7712;-diag-disable 5194>")
if (${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel*")
  set(CMAKE_Fortran_FLAGS_DEBUG  "${CMAKE_Fortran_FLAGS_DEBUG} -g -traceback -check")
  add_compile_options("-fPIE")
  add_compile_options("-fPIC")
  # -recusrive is needed to ensure that routines for accessing PDFs and
  # coupling are thread-safe
  # (also tried  -reentrancy:threaded, but it didn't work
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -recursive")
endif()


# NVidia
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,NVHPC>:-Minform=warn;>")
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,NVHPC>:-Mnofpapprox>")
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,NVHPC>:-target=multicore;-stdpar=multicore>")
if (${CMAKE_Fortran_COMPILER_ID} MATCHES "NVHPC")
  #message(WARNING "HOPPET: ${CMAKE_Fortran_COMPILER_ID} ${COMPILER_BUG_MSG}")
  set(CMAKE_Fortran_FLAGS_DEBUG  "${CMAKE_Fortran_FLAGS_DEBUG} -g -traceback -Mbounds -Mchkptr")
endif()

# NAG
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,NAG>:-dcfuns;-dusty>")
if (${CMAKE_Fortran_COMPILER_ID} MATCHES "NAG")
  set(CMAKE_Fortran_FLAGS_DEBUG  "${CMAKE_Fortran_FLAGS_DEBUG} -g -C=all -mtrace -gline")
endif()

# IBM XL
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,XL>:-qfree;-qextname;-qhalt=s;-qxlf2003=polymorphic>")
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,XL>:-qstackprotect=all>")
if (${CMAKE_Fortran_COMPILER_ID} MATCHES "XL")
  # xlf seems to be working OK in CI, so comment out the COMPILER_BUG_MSG
  #message(WARNING "HOPPET: ${CMAKE_Fortran_COMPILER_ID} ${COMPILER_BUG_MSG}")
  # GPS 2025-02-22: added -qxlf2003=polymorphic to be able to use F2003 features needed in mass_thresholds_n3lo.f90
  # but this does not seem to solve the problem
  set(CMAKE_Fortran_FLAGS_DEBUG  "${CMAKE_Fortran_FLAGS_DEBUG} -g -qcheck -qsigtrap")
endif()

# AMD Flang
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,Flang>:-Wall>")
if (${CMAKE_Fortran_COMPILER_ID} MATCHES "Flang")
  #message(WARNING "HOPPET: ${CMAKE_Fortran_COMPILER_ID} ${COMPILER_BUG_MSG}")
  set(CMAKE_Fortran_FLAGS_DEBUG  "${CMAKE_Fortran_FLAGS_DEBUG} -g")
endif()

# ARM 
if (${CMAKE_Fortran_COMPILER_ID} MATCHES "ARM")
  #message(WARNING "HOPPET: ${CMAKE_Fortran_COMPILER_ID} ${COMPILER_BUG_MSG}")
  set(CMAKE_Fortran_FLAGS_DEBUG  "${CMAKE_Fortran_FLAGS_DEBUG} -g")
endif()

if (HOPPET_ENABLE_FPES)
  #add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-ffpe-trap=zero,overflow,underflow,invalid,denormal>")
  # note: underflow and denormal are potentially problematic with some parts of the code
  # e.g. if they do an expansion in a high-order polynomial, the result may be zero
  add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-ffpe-trap=zero,overflow,invalid>")
  add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-fpe0>")
  add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>:-fpe0>")
  add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>:-O3>") # Intel needs explicit -O3 with the debug flag
  add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,NVHPC>:-Ktrap=divz,denorm,fp,inexact,inv,ovf,unf>")
  add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,NAG>:-ieee=stop>")
  add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,XL>:-qflttrap=overflow:underflow:invalid:zerodivide:enable>")
  add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,Flang>:-ffp-exception-behavior=strict>")
  add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,ARM>:-ffp-contract=on>")
endif()

# test for the array bound bug
if (${CMAKE_Fortran_COMPILER_ID} MATCHES "XL")
  message(WARNING "HOPPET: XLF-compiled executable may hang on array bound test, skipping; make sure you run ctest")
else()
  SET(ARRAY_BOUND_TEST "${PROJECT_SOURCE_DIR}/cmake/Tests/array_bound_bug.f90")
  try_run(ARRAY_BOUND_BUG_RUN ARRAY_BOUND_BUG_COMPILE
    ${CMAKE_BINARY_DIR} ${ARRAY_BOUND_TEST}
    COMPILE_OUTPUT_VARIABLE COMPILE_OUTPUT
    RUN_OUTPUT_VARIABLE RUN_OUTPUT
    )
  if (NOT ARRAY_BOUND_BUG_COMPILE)
    # options includes SEND_ERROR or WARNING
    message(WARNING "HOPPET: Failed to compile ${ARRAY_BOUND_TEST}")
    message(WARNING "Output was ${COMPILE_OUTPUT}")
  elseif(NOT ARRAY_BOUND_BUG_RUN EQUAL 0)
    set (COMPILER_BUG_MSG "Compiler appears to contain bugs that prevent HOPPET from working correctly. Please see https://fortran-lang.discourse.group/t/an-interesting-difference-between-compilers/7131 for details and one of the bug reports at https://forums.developer.nvidia.com/t/bugreport-for-nvfortran-treatment-of-array-bounds-in-the-specification-part/278324")
    message(WARNING "HOPPET: ${CMAKE_Fortran_COMPILER_ID} compilation & run of ${ARRAY_BOUND_TEST} gave incorrect results")
    message(WARNING "Output was ${RUN_OUTPUT}")
    message(WARNING "for more info: ${COMPILER_BUG_MSG}")
  else()
    message(STATUS "HOPPET: ${CMAKE_Fortran_COMPILER_ID} compiler passes array bound test")
  endif()
endif()

set(HOPPET_BUILD_PYINTERFACE_NOCACHE ${HOPPET_BUILD_PYINTERFACE})
if (HOPPET_BUILD_PYINTERFACE OR HOPPET_BUILD_PYWHEELS)
  find_package(SWIG)
  if (SWIG_FOUND)
    message(STATUS "HOPPET: SWIG found: ${SWIG_EXECUTABLE}")
  else()
    message(WARNING "HOPPET: Python interface requested, but SWIG not found, disabling Python interface")
    set(HOPPET_BUILD_PYINTERFACE_NOCACHE OFF) 
  endif()
  message(STATUS "HOPPET: SWIG_VERSION=${SWIG_VERSION} SWIG_LIBRARIES=${SWIG_LIBRARIES} SWIG_INCLUDE_DIRS=${SWIG_INCLUDE_DIRS}")
endif()


# Output of config info   ==============================================
message(STATUS "HOPPET: CMAKE_Fortran_COMPILER_ID=${CMAKE_Fortran_COMPILER_ID}")
message(STATUS "HOPPET: CMAKE_Fortran_FLAGS, CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE} ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}}" )
message(STATUS "HOPPET: CMAKE_C_COMPILER_ID=${CMAKE_C_COMPILER_ID}")
message(STATUS "HOPPET: CMAKE_CXX_COMPILER_ID=${CMAKE_CXX_COMPILER_ID}")
message(STATUS "HOPPET: HOPPET_USE_EXACT_COEF=${HOPPET_USE_EXACT_COEF}")
message(STATUS "HOPPET: HOPPET_BUILD_EXAMPLES=${HOPPET_BUILD_EXAMPLES}")
message(STATUS "HOPPET: HOPPET_BUILD_PYINTERFACE=${HOPPET_BUILD_PYINTERFACE_NOCACHE}")
message(STATUS "HOPPET: HOPPET_BUILD_PYWHEELS=${HOPPET_BUILD_PYWHEELS}")
message(STATUS "HOPPET: HOPPET_ENABLE_TESTING=${HOPPET_ENABLE_TESTING}")
message(STATUS "HOPPET: HOPPET_BUILD_BENCHMARK=${HOPPET_BUILD_BENCHMARK}")
message(STATUS "HOPPET: HOPPET_BUILD_TESTDIR=${HOPPET_BUILD_TESTDIR}")
message(STATUS "HOPPET: HOPPET_ENABLE_FPES=${HOPPET_ENABLE_FPES}")
message(STATUS "HOPPET: HOPPET_ENABLE_DEBUG: ${HOPPET_ENABLE_DEBUG}")


if (HOPPET_BUILD_EXAMPLES)
  find_package(LHAPDF) 
  if (LHAPDF_FOUND)
    message(STATUS "HOPPET: LHAPDF_VERSION=${LHAPDF_VERSION} LHAPDF_LIBRARIES=${LHAPDF_LIBRARIES} LHAPDF_INCLUDE_DIRS=${LHAPDF_INCLUDE_DIRS}")
  else()
    message(WARNING "HOPPET: LHAPDF not found, examples that need it will not be built")
  endif()
endif()


# Sources ==============================================================
file(GLOB SOURCESf90 src/*.f90 src/*.F90 src/splitting-functions/*.f90 src/mass-thresholds/*.f90)
file(GLOB SOURCESf src/*.f src/splitting-functions/*.f src/mass-thresholds/*.f src/param-coefs/*.f)

#------------------------------------------------------------------------- 
# Define the two required variables before including
# the source code for watching a git repository.
option(GIT_WATCHER "Add source code which allows runs for watching the status of the git repository" ON)
set(PRE_CONFIGURE_FILE "src/config/hoppet_git_state.f90.in")
set(POST_CONFIGURE_FILE "${CMAKE_CURRENT_BINARY_DIR}/hoppet_git_state.f90")
if (GIT_WATCHER)
  message("-- Enabling the git watcher support")
  set(GIT_FAIL_IF_NONZERO_EXIT FALSE)
  include(cmake/git_watcher.cmake)  
else()
  message("-- Disabling the git watcher support")  
  # just copy PRE_CONFIGURE_FILE to POST_CONFIGURE_FILE
  configure_file(${PRE_CONFIGURE_FILE} ${POST_CONFIGURE_FILE} COPYONLY)
endif()
list(APPEND SOURCESf90 ${POST_CONFIGURE_FILE})


# alternative is a manual list (NB: this specific list is not complete)
# set(SOURCES
#   hplog.f xpij2e.f xpns2e.f assertions.f90 coefficient_functions.f90
#   convolution.f90 dglap_choices.f90 dglap_holders.f90 dglap_objects.f90
#   dummy_pdfs.f90 evolution.f90 hoppet.f90 integrator.f90
#   interpolation.f90 new_as.f90 pdf_general.f90 pdf_representation.f90
#   pdf_tabulate.f90 qcd.f90 qcd_coupling.f90 random.f90 runge_kutta.f90
#   sort.f90 special_functions.f90 splitting_functions.f90
#   splitting_functions_nnlo.f90 streamlined_interface.f90 types.f90
#   warnings_and_errors.f90 welcome_message.f90 xa2hgp.f90 xpij2n.f90
#   xpij2p.f90 xpns2n.f90 xpns2p.f90)

# and additionally in the directories used for conditional compilation
file(GLOB EXACT_COEF_SOURCESf90 src/exact-coefs/*.f90)
file(GLOB EXACT_COEF_SOURCESf src/exact-coefs/*.f)
file(GLOB NO_EXACT_COEF_SOURCES src/no-exact-coefs/*.f src/no-exact-coefs/*.f90)
file(GLOB NO_EXACT_COEF_SOURCESf90 src/no-exact-coefs/*.f90)
file(GLOB NO_EXACT_COEF_SOURCESf src/no-exact-coefs/*.f)


if (${CMAKE_Fortran_COMPILER_ID} MATCHES "XL") 
 set_source_files_properties( ${SOURCESf}  PROPERTIES COMPILE_OPTIONS "-qfixed")
endif()

# if HOPPET_USE_EXACT_COEF is ON, use the exact coefficient functions sources
if(HOPPET_USE_EXACT_COEF)
  set(SOURCESf90 ${SOURCESf90} ${EXACT_COEF_SOURCESf90})
  set(SOURCESf ${SOURCESf} ${EXACT_COEF_SOURCESf})
else()
  set(SOURCESf90 ${SOURCESf90} ${NO_EXACT_COEF_SOURCESf90})
  set(SOURCESf ${SOURCESf} ${NO_EXACT_COEF_SOURCESf})
endif()


#  Libraries ===========================================================
add_library(hoppet SHARED ${SOURCESf90} ${SOURCESf})
target_include_directories(hoppet PUBLIC $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/hoppet> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hoppet> )
set_target_properties(hoppet PROPERTIES OUTPUT_NAME hoppet
    ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/$<0:>
    LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/$<0:>
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/$<0:>
    Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/hoppet_shared
    SOVERSION 0.0.0
    EXPORT_NAME hoppet::hoppet
    )
add_library(hoppet::hoppet ALIAS hoppet)


add_library(hoppet_static STATIC ${SOURCESf90} ${SOURCESf})
target_include_directories(hoppet_static PUBLIC $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/hoppet> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hoppet> )
set_target_properties(hoppet_static PROPERTIES OUTPUT_NAME hoppet
    ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/$<0:>
    LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/$<0:>
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/$<0:>
    Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/hoppet
    EXPORT_NAME hoppet::hoppet_static
    )
add_library(hoppet::hoppet_static ALIAS hoppet_static)

if (HOPPET_ENABLE_TESTING)
  ENABLE_TESTING()
  add_subdirectory(unit-tests)
endif()

# Examples =============================================================
#
# Note that below, all examples link against the static library.
# On Linux this avoids the risk of accidentally loading an installed
# version of the library when the user has an LD_LIBRARY_PATH set, 
# because LD_LIBRARY_PATH always takes precedence over the RUNPATH
# that CMake sets in the executable.
#
# Building the static library does not seem to add much time to the
# build (a few percent)
if (HOPPET_BUILD_EXAMPLES)
  set( f90programs_f90 
                   structure_functions_example  
                   structure_functions_example_flavour  
                   tabulation_example_2                  
                   tabulation_example      
                   tabulation_example_n3lo
                   tabulation_example_qed_streamlined
                   tabulation_example_up_and_down
                   tabulation_example_diehl_crosschecks  
                   tabulation_example_qed  
                   tabulation_example_streamlined
                   structure_functions_example_n3lo_tests
		   )
  # first handle the f90 examples
  foreach ( p ${f90programs_f90})
    add_executable( ${p}_f90 example_f90/${p}.f90)
    target_link_libraries(${p}_f90 PUBLIC hoppet::hoppet_static)
    set_target_properties(${p}_f90 PROPERTIES OUTPUT_NAME ${p}
      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/example_f90/
      Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/example_f90/
      POSITION_INDEPENDENT_CODE TRUE
    )
  endforeach()
  if (LHAPDF_FOUND)
    set( f90programs_lhapdf fast_pdf_evaluation fast_manypdf_evaluation)
    foreach ( p ${f90programs_lhapdf})
      add_executable( ${p}_f90 example_f90/with-lhapdf/${p}.f90 example_f90/with-lhapdf/hoppet_lhapdf.f90)
      target_link_libraries(${p}_f90 PUBLIC hoppet::hoppet_static PUBLIC LHAPDF::LHAPDF)
      set_target_properties(${p}_f90 PROPERTIES OUTPUT_NAME ${p}
        RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/example_f90/
        Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/example_f90/modules/${p}_modules
        POSITION_INDEPENDENT_CODE TRUE
        )
      if (${CMAKE_Fortran_COMPILER_ID} MATCHES "XL") 
        set_source_files_properties( example_f90/${p}.f  PROPERTIES COMPILE_OPTIONS "-qfixed")
      endif()     
    endforeach()
  endif()


  # then the f77 examples
  set( f77programs_f tabulation_example)
  if (LHAPDF_FOUND)
    list(APPEND f77programs_f compare_lhapdf_hoppet convolution_example)  
  endif()
  foreach ( p ${f77programs_f})
    add_executable( ${p}_f77 example_f77/${p}.f)
    target_link_libraries(${p}_f77 PUBLIC hoppet::hoppet_static)
    set_target_properties(${p}_f77 PROPERTIES  OUTPUT_NAME ${p} 
      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/example_f77/
      Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/example_f77/
      POSITION_INDEPENDENT_CODE TRUE
      )
    if (${CMAKE_Fortran_COMPILER_ID} MATCHES "XL") 
      set_source_files_properties( example_f77/${p}.f  PROPERTIES COMPILE_OPTIONS "-qfixed")
    endif()     
  endforeach()
  if (LHAPDF_FOUND)
    target_link_libraries(compare_lhapdf_hoppet_f77 PUBLIC LHAPDF::LHAPDF)
    target_link_libraries(convolution_example_f77 PUBLIC LHAPDF::LHAPDF)
  endif()

  # and finally the C++ examples
  set( cppprograms structure_functions_example tabulation_example_qed tabulation_example)
  foreach ( p ${cppprograms})
    add_executable( ${p}_cpp example_cpp/${p}.cc)
    target_include_directories(${p}_cpp PUBLIC ${PROJECT_SOURCE_DIR}/src)
    target_link_libraries(${p}_cpp PUBLIC hoppet::hoppet_static)
    set_target_properties(${p}_cpp PROPERTIES  OUTPUT_NAME ${p}
      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/example_cpp/
      Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/example_cpp/
      POSITION_INDEPENDENT_CODE TRUE
     )
  endforeach()
  if (LHAPDF_FOUND)
    set(cppprograms_lhapdf fast_pdf_evaluation)
    foreach ( p ${cppprograms_lhapdf})
      add_executable( ${p}_cpp example_cpp/with-lhapdf/${p}.cc)
      target_include_directories(${p}_cpp PUBLIC ${PROJECT_SOURCE_DIR}/src)
      target_link_libraries(${p}_cpp PUBLIC hoppet::hoppet_static PUBLIC LHAPDF::LHAPDF)
      set_target_properties(${p}_cpp PROPERTIES  OUTPUT_NAME ${p}
                            RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/example_cpp/
                            Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/example_cpp/
                            POSITION_INDEPENDENT_CODE TRUE
      )
    endforeach()
  endif()  

  if (HOPPET_ENABLE_TESTING)
    set( execs 
      tabulation_example 
      tabulation_example_n3lo 
      tabulation_example_qed 
      structure_functions_example
      structure_functions_example_flavour
      tabulation_example_qed_streamlined
      #tabulation_example_up_and_down 
      )
    if(HOPPET_USE_EXACT_COEF)
      list(APPEND execs structure_functions_example:exact )
    endif()
    foreach (ex ${execs})
      set(testName ${ex})
      if (ex MATCHES "tabulation_example_n3lo")
        set(testName "${testName}_takes_up_to_30s")
      endif()
      if (ex MATCHES "structure_functions_example:exact")
        set(testName "${testName}_takes_up_to_30s")
      endif()
      add_test(NAME ${testName} COMMAND 
         sh -c "${CMAKE_SOURCE_DIR}/scripts/check example_f90/${ex}")
    endforeach()
    set( execs_cpp tabulation_example)
    foreach (ex ${execs_cpp})
      set(testName "${ex}_cpp")
      add_test(NAME ${testName} COMMAND 
         sh -c "${CMAKE_SOURCE_DIR}/scripts/check example_cpp/${ex}")
    endforeach()
  endif() # end of HOPPET_ENABLE_TESTING

endif()

if (HOPPET_BUILD_BENCHMARK)
  set( bms benchmark prec_and_timing structure_functions_benchmark_checks 
           evolution_n3lo_tests sumrules compare2files_v2) 
  foreach ( b ${bms})
    add_executable( ${b} benchmarking/${b}.f90 benchmarking/io_utils.f90 benchmarking/lcl_dec.f90 benchmarking/NameSelect.f90)
    target_link_libraries(${b} PUBLIC hoppet::hoppet_static)
    #target_link_libraries(${b} PUBLIC LHAPDF::LHAPDF)
    set_target_properties(${b} PROPERTIES  OUTPUT_NAME ${b}
      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/benchmarking/
      # each target gets its modules built separately, even though some 
      # could in principle be shared (which would require a different CMake setup)
      Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/benchmarking/modules/${b}_modules
      POSITION_INDEPENDENT_CODE TRUE
     )
  endforeach()
  if (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
    # turn off a spurious warning in benchmarking/benchmark.f90
    target_compile_options(benchmark PRIVATE "-Wno-do-subscript")
  endif()

  if (HOPPET_ENABLE_TESTING)
    # this list of tests should align with the list of runs at the start of
    # benchmarking/prec_and_timing.default_output/gen-default-output.sh
    set (prec_and_timing_tests 
      benchmarking/prec_and_timing:-dy:0.2:-olnlnQ:2:-yinterp-order:2
      benchmarking/prec_and_timing:-dy:0.2:-olnlnQ:3:-yinterp-order:3
      benchmarking/prec_and_timing:-dy:0.2:-olnlnQ:4:-yinterp-order:4
      benchmarking/prec_and_timing:-dy:0.1:-exact-nnlo-th:-exact-nnlo-sp
    )
    foreach (test ${prec_and_timing_tests})
      set(testName ${test})
      add_test(NAME ${testName} COMMAND 
         sh -c "${CMAKE_SOURCE_DIR}/scripts/check ${test}")
    endforeach()
  endif()

  if (LHAPDF_FOUND)
    set(lhapdf_cpp_bms lhapdf_timings)
    foreach ( p ${lhapdf_cpp_bms})
      add_executable( ${p}_cpp benchmarking/${p}.cc)
      target_link_libraries(${p}_cpp PUBLIC hoppet::hoppet_static PUBLIC LHAPDF::LHAPDF)
      set_target_properties(${p}_cpp PROPERTIES  OUTPUT_NAME ${p}
                            RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/benchmarking
                            POSITION_INDEPENDENT_CODE TRUE
      )
    endforeach()
  endif()  


endif()

if (HOPPET_BUILD_TESTDIR)
  set( testcode bidirectional_mass_thresholds coupling-fixed coupling determine_accuracies intrinsic_charm
    lumi_tests msbar_tests split_mat_ops streamlined_multisplit test_qed test_qed_lepton_evol test_qed_obj
    timelike_NLOevol_tests timelike_tests test_split_integration test-n3lo-mtm) 
  foreach ( b ${testcode})
    add_executable( ${b} tests/${b}.f90 benchmarking/io_utils.f90 benchmarking/lcl_dec.f90 benchmarking/NameSelect.f90)
    target_link_libraries(${b} PUBLIC hoppet::hoppet_static)
    set_target_properties(${b} PROPERTIES  OUTPUT_NAME ${b}
      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/tests/
      Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/tests/${b}_modules
      POSITION_INDEPENDENT_CODE TRUE
      )
  endforeach()
  if(LHAPDF_FOUND)
    set( testcode_lhapdf tabulation_n3lo_v130_paper  test_qed_evol_lhapdf  test_qed_lhapdf)
    foreach ( b ${testcode_lhapdf})
      add_executable( ${b} tests/with-lhapdf/${b}.f90 benchmarking/io_utils.f90 benchmarking/lcl_dec.f90 benchmarking/NameSelect.f90)
      target_link_libraries(${b} PUBLIC hoppet::hoppet_static)
      target_link_libraries(${b} PUBLIC LHAPDF::LHAPDF)
      set_target_properties(${b} PROPERTIES  OUTPUT_NAME ${b}
	RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/tests/
	Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/tests/${b}_modules
	POSITION_INDEPENDENT_CODE TRUE
	)
    endforeach()
  endif()
endif()

if (HOPPET_BUILD_PYINTERFACE_NOCACHE OR HOPPET_BUILD_PYWHEELS)
  add_subdirectory(pyinterface)
endif()

# Cmake config files   #################################################
include(CMakePackageConfigHelpers)

set(CONFIG_INSTALL_DIR ${CMAKE_INSTALL_DATADIR}/hoppet/cmake)

configure_package_config_file(cmake/Templates/hoppetConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/share/hoppet/cmake/hoppetConfig.cmake
  INSTALL_DESTINATION ${CONFIG_INSTALL_DIR}
  PATH_VARS CMAKE_INSTALL_INCLUDEDIR CMAKE_INSTALL_LIBDIR)

write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/share/hoppet/cmake/hoppetConfig-version.cmake COMPATIBILITY SameMajorVersion VERSION ${PROJECT_VERSION})

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/share/hoppet/cmake/hoppetConfig.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/share/hoppet/cmake/hoppetConfig-version.cmake
              DESTINATION ${CONFIG_INSTALL_DIR} COMPONENT devel)

# Install  #############################################################
# Create symlink for backwards compatibility
#install(CODE "
#  set(libdir \"\$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}\")
#  file(CREATE_LINK \"libhoppet.a\" \"\${libdir}/libhoppet_v1.a\" SYMBOLIC)
#")
install(TARGETS hoppet hoppet_static  EXPORT hoppetTargets DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/hoppet DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ )
install(FILES ${PROJECT_SOURCE_DIR}/src/hoppet.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ )
install(FILES ${PROJECT_BINARY_DIR}/bin/hoppet-config 
                                DESTINATION ${CMAKE_INSTALL_BINDIR} 
                                PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

install(EXPORT hoppetTargets DESTINATION ${CONFIG_INSTALL_DIR} COMPONENT devel)

# hoppet-config script =================================================
set(prefix ${CMAKE_INSTALL_PREFIX})
set(modprefix "${prefix}/${CMAKE_INSTALL_INCLUDEDIR}/hoppet")
set(VERSION "${HOPPET_VERSION}")
set(libdir "${CMAKE_INSTALL_LIBDIR}")
configure_file("${PROJECT_SOURCE_DIR}/scripts/hoppet-config.in" "${PROJECT_BINARY_DIR}/bin/hoppet-config")

# Banner to print after install
install(CODE "execute_process(COMMAND bash \"${CMAKE_SOURCE_DIR}/scripts/print_banner.sh\" ${PROJECT_VERSION})")

