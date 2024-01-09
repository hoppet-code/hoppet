
## mkdir build
## cmake -S . -B build [-DCMAKE_INSTALL_PREFIX=/opt/myproject] [-DHOPPET_USE_EXACT_COEF=ON] ..
## cmake --build  build -j 
## cmake --install build
## ctest --test-dir build

cmake_minimum_required(VERSION 3.15)
# set the project name
project(hoppet VERSION 1.3.0 LANGUAGES C CXX Fortran)
set(hoppet_VERSION 1.3.0-devel)
set(PROJECT_VERSION 1.3.0-devel)
set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules ${CMAKE_MODULE_PATH})
include("GNUInstallDirs")

# Add options ==========================================================
option(HOPPET_USE_EXACT_COEF  "Use exact coefficient functions"  OFF)
option(HOPPET_BUILD_EXAMPLES  "Build examples" ON)
option(HOPPET_ENABLE_TESTING  "Enable testing. Requires building the examples." ON)
option(HOPPET_BUILD_BENCHMARK "Build benchmark." ON)
option(HOPPET_ENABLE_FPES     "Enable trapping for the floating point exceptions." OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RELEASE")
  message(STATUS "HOPPET: Build type not specified, use: ${CMAKE_BUILD_TYPE}")
endif()

string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE)
message(STATUS "HOPPET: CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")

# Flags for compilers ==============================================
# GNU
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-Wall>")
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-Wextra;-Wno-unused-parameter;-Wno-unused-variable>")
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-Wno-unused-dummy-argument;-Wno-maybe-uninitialized>")
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-Wno-unused;-Wno-uninitialized;-Wno-conversion>")
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-Wno-compare-reals;-Wno-function-elimination>")
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-Wno-integer-division>")
if (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS_DEBUG  "${CMAKE_Fortran_FLAGS_DEBUG} -g -fbacktrace -fcheck=all")
endif()

# Intel
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-diag-disable 8291;-diag-disable 10448;-diag-disable 8889;-diag-disable 6717;-diag-disable 8291;-diag-disable 7712;-diag-disable 5194>")
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>:-warn;-diag-disable 8291;-diag-disable 8889;-diag-disable 6717;-diag-disable 8291;-diag-disable 7712;-diag-disable 5194>")
if (${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel*")
  set(CMAKE_Fortran_FLAGS_DEBUG  "${CMAKE_Fortran_FLAGS_DEBUG} -g -traceback -check")
endif()

# NVidia
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,NVHPC>:-Minform=warn;>")
if (${CMAKE_Fortran_COMPILER_ID} MATCHES "NVHPC")
  set(CMAKE_Fortran_FLAGS_DEBUG  "${CMAKE_Fortran_FLAGS_DEBUG} -g -traceback -Mbounds -Mchkptr")
endif()

# NAG
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,NAG>:-dcfuns;-dusty>")
if (${CMAKE_Fortran_COMPILER_ID} MATCHES "NAG")
  set(CMAKE_Fortran_FLAGS_DEBUG  "${CMAKE_Fortran_FLAGS_DEBUG} -g -C=all -mtrace -gline")
endif()

# XL
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,XL>:-qfree;-qextname;-qhalt=s>")
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,XL>:-qstackprotect=all>")

if (${CMAKE_Fortran_COMPILER_ID} MATCHES "XL")
  set(CMAKE_Fortran_FLAGS_DEBUG  "${CMAKE_Fortran_FLAGS_DEBUG} -g -qcheck -qsigtrap")
endif()

if (HOPPET_ENABLE_FPES)
  add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-ffpe-trap=zero,overflow,underflow,invalid,denormal>")
  add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-fpe0>")
  add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>:-fpe0>")
  add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,NVHPC>:-Ktrap=divz,denorm,fp,inexact,inv,ovf,unf>")
  add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,NAG>:-ieee=stop>")
  add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,XL>:-qflttrap=overflow:underflow:invalid:zerodivide:enable>")
endif()


# Output of config info   ==============================================
message(STATUS "HOPPET: CMAKE_Fortran_COMPILER_ID=${CMAKE_Fortran_COMPILER_ID}")
message(STATUS "HOPPET: CMAKE_Fortran_FLAGS, CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE} ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}}" )
message(STATUS "HOPPET: CMAKE_C_COMPILER_ID=${CMAKE_C_COMPILER_ID}")
message(STATUS "HOPPET: CMAKE_CXX_COMPILER_ID=${CMAKE_CXX_COMPILER_ID}")
message(STATUS "HOPPET: HOPPET_USE_EXACT_COEF=${HOPPET_USE_EXACT_COEF}")
message(STATUS "HOPPET: HOPPET_BUILD_EXAMPLES=${HOPPET_BUILD_EXAMPLES}")
message(STATUS "HOPPET: HOPPET_ENABLE_TESTING=${HOPPET_ENABLE_TESTING}")
message(STATUS "HOPPET: HOPPET_BUILD_BENCHMARK=${HOPPET_BUILD_BENCHMARK}")
message(STATUS "HOPPET: HOPPET_ENABLE_FPES=${HOPPET_ENABLE_FPES}")
find_package(LHAPDF) 
message(STATUS "HOPPET: LHAPDF_VERSION=${LHAPDF_VERSION} LHAPDF_LIBRARIES=${LHAPDF_LIBRARIES} LHAPDF_INCLUDE_DIRS=${LHAPDF_INCLUDE_DIRS}")


# hoppet-config script =================================================
set(prefix ${CMAKE_INSTALL_PREFIX})
set(modprefix "${prefix}/${CMAKE_INSTALL_INCLUDEDIR}/hoppet")
set(VERSION "${hoppet_VERSION}")
configure_file("${PROJECT_SOURCE_DIR}/hoppet-config.in" "${PROJECT_BINARY_DIR}/bin/hoppet-config")


# Sources ==============================================================
file(GLOB SOURCESf90 src/*.f90)
file(GLOB SOURCESf src/*.f)


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
file(GLOB EXACT_COEF_SOURCESf90 src/exact-coefs/*.f90)
file(GLOB EXACT_COEF_SOURCESf src/exact-coefs/*.f)
file(GLOB NO_EXACT_COEF_SOURCES src/no-exact-coefs/*.f src/no-exact-coefs/*.f90)
file(GLOB NO_EXACT_COEF_SOURCESf90 src/no-exact-coefs/*.f90)
file(GLOB NO_EXACT_COEF_SOURCESf src/no-exact-coefs/*.f90)


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
set_target_properties(hoppet PROPERTIES OUTPUT_NAME hoppet_v1
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
set_target_properties(hoppet_static PROPERTIES OUTPUT_NAME hoppet_v1
    ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/$<0:>
    LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/$<0:>
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/$<0:>
    Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/hoppet
    EXPORT_NAME hoppet::hoppet_static
    )
add_library(hoppet::hoppet_static ALIAS hoppet_static)

# Examples =============================================================
if (HOPPET_BUILD_EXAMPLES)
  set( f90programs_f90 
                   structure_functions_example  
                   tabulation_example_2                  
                   tabulation_example      
                   tabulation_example_qed_streamlined
                   sumrules                     
                   tabulation_example_diehl_crosschecks  
                   tabulation_example_qed  
                   tabulation_example_streamlined
       )
  set( f90programs_f hplog5_tests )
  foreach ( p ${f90programs_f90})
    add_executable( ${p}_f90 example_f90/${p}.f90)
    target_link_libraries(${p}_f90 PUBLIC hoppet::hoppet_static)
    set_target_properties(${p}_f90 PROPERTIES OUTPUT_NAME ${p}
      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/example_f90/
      Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/example_f90/
      POSITION_INDEPENDENT_CODE TRUE
     )
  endforeach()
  foreach ( p ${f90programs_f})
    add_executable( ${p}_f90 example_f90/${p}.f)
    target_link_libraries(${p}_f90 PUBLIC hoppet::hoppet_static)
    set_target_properties(${p}_f90 PROPERTIES OUTPUT_NAME ${p}
      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/example_f90/
      Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/example_f90/
      POSITION_INDEPENDENT_CODE TRUE
     )
    if (${CMAKE_Fortran_COMPILER_ID} MATCHES "XL") 
      set_source_files_properties( example_f90/${p}.f  PROPERTIES COMPILE_OPTIONS "-qfixed")
    endif()     
  endforeach()
  set( f77programs_f
       tabulation_example
       xc2ns3e_vs_xc2ns3p)
 if (LHAPDF_FOUND)
   list(APPEND f77programs_f compare_lhapdf_hoppet convolution_example)  
 endif()
 set( f77programs_c cpp_structure_functions_example  cpp_tabulation_example_qed cpp_tabulation_example )
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
  foreach ( p ${f77programs_c})
    add_executable( ${p}_f77c example_f77/${p}.cc)
    target_link_libraries(${p}_f77c PUBLIC hoppet::hoppet_static)
    set_target_properties(${p}_f77c PROPERTIES  OUTPUT_NAME ${p}
      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/example_f77/
      Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/example_f77/
      POSITION_INDEPENDENT_CODE TRUE
     )
  endforeach()
  if (HOPPET_ENABLE_TESTING)
    ENABLE_TESTING()
    set( execs tabulation_example tabulation_example_qed structure_functions_example tabulation_example_qed_streamlined)
    foreach ( ex ${execs})
      add_test(NAME ${ex}_copy  COMMAND sh -c "cat ${PROJECT_SOURCE_DIR}/example_f90/${ex}.default_output | grep -vi -e '^ *[a-z]' -e '^ *$' > ${PROJECT_BINARY_DIR}/example_f90/${ex}.default_output")
      add_test(NAME ${ex} COMMAND sh -c "${PROJECT_BINARY_DIR}/example_f90/${ex} | grep -vi -e '^ *[a-z]' -e '^ *$' > ${PROJECT_BINARY_DIR}/example_f90/${ex}.this_output" )
      add_test(NAME ${ex}_compare  COMMAND ${CMAKE_COMMAND} -E compare_files ${PROJECT_BINARY_DIR}/example_f90/${ex}.this_output  ${PROJECT_BINARY_DIR}/example_f90/${ex}.default_output)
      set_tests_properties( ${ex}_compare PROPERTIES DEPENDS "${ex};${ex}_copy")
    endforeach()    
  endif()
endif()
if (HOPPET_BUILD_BENCHMARK)
  set( bms benchmark prec_and_timing structure_functions_benchmark_checks) 
  foreach ( b ${bms})
    add_executable( ${b} benchmarking/${b}.f90 benchmarking/io_utils.f90 benchmarking/lcl_dec.f90 benchmarking/NameSelect.f90)
    target_link_libraries(${b} PUBLIC hoppet::hoppet_static)
    set_target_properties(${b} PROPERTIES  OUTPUT_NAME ${b}
      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/benchmarking/
      Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/benchmarking/${b}_modules
      POSITION_INDEPENDENT_CODE TRUE
     )
  endforeach()

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
install(TARGETS hoppet hoppet_static  EXPORT hoppetTargets DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/hoppet DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ )
install(FILES ${PROJECT_SOURCE_DIR}/src/hoppet_v1.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ )
install(FILES ${PROJECT_BINARY_DIR}/bin/hoppet-config 
                                DESTINATION ${CMAKE_INSTALL_BINDIR} 
                                PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
install(EXPORT hoppetTargets DESTINATION ${CONFIG_INSTALL_DIR} COMPONENT devel)
