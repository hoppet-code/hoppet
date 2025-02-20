cmake_minimum_required(VERSION 3.15)
project(pyinterface LANGUAGES C CXX)

#include(FortranCInterface)
# these are useful for more advanced work?
#FortranCInterface_VERIFY()
#FortranCInterface_HEADER(${CMAKE_BINARY_DIR}/FortranCInterface.h)


# Find Python and SWIG
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

# # Find Hoppet
# find_file(HOPPET_CONFIG hoppet-config)
# if(NOT HOPPET_CONFIG)
#     message(FATAL_ERROR "Could not find hoppet-config")
# endif()
# 
# execute_process(
#   COMMAND bash -c "${HOPPET_CONFIG} --prefix"
#   OUTPUT_VARIABLE HOPPET_PREFIX
#   OUTPUT_STRIP_TRAILING_WHITESPACE
# )
# 
# execute_process(
#   COMMAND bash -c "${HOPPET_CONFIG} --libs"
#   OUTPUT_VARIABLE HOPPET_LIBRARY
#   OUTPUT_STRIP_TRAILING_WHITESPACE
# )
# 
# set(HOPPET_INCLUDE_DIR "${HOPPET_PREFIX}/include/hoppet")
# set(HOPPET_LIB_DIR "${HOPPET_PREFIX}/lib")
# set(HOPPET_LIBRARY "${HOPPET_LIBRARY} -Wl,-rpath,${HOPPET_LIB_DIR}")

# Include directories
include_directories(${Python3_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/src)

# SWIG interface
set(SWIG_INPUT_FILE hoppet_v1.i)
set(CMAKE_SWIG_FLAGS -doxygen)

set_source_files_properties(${SWIG_INPUT_FILE} PROPERTIES CPLUSPLUS ON)

# Generate the SWIG wrapper
swig_add_library(hoppet_v1
  TYPE MODULE
  LANGUAGE python
  SOURCES ${SWIG_INPUT_FILE}
)

# Ensure the SWIG-generated wrapper file is compiled as C++
set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/hoppet_v1PYTHON_wrap.cxx PROPERTIES LANGUAGE CXX)

# Set C++ standard
set_target_properties(${SWIG_MODULE_hoppet_v1_REAL_NAME} PROPERTIES CXX_STANDARD 11)

# Link libraries
target_link_libraries(${SWIG_MODULE_hoppet_v1_REAL_NAME}
#PRIVATE ${HOPPET_LIBRARY} ${Fortran_LIBRARIES} ${Python3_LIBRARIES}
  PUBLIC hoppet::hoppet ${Fortran_LIBRARIES} ${Python3_LIBRARIES}
#PUBLIC hoppet_v1 ${Fortran_LIBRARIES} ${Python3_LIBRARIES}
)

# Install the Python module
#message("******* ${SWIG_MODULE_hoppet_v1_REAL_NAME}")
install(TARGETS ${SWIG_MODULE_hoppet_v1_REAL_NAME}
  LIBRARY DESTINATION ${Python3_SITEARCH}
)
# manually install hoppet_v1.py -- not sure why this is needed
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/hoppet_v1.py
  DESTINATION ${Python3_SITEARCH}
)

# Add the tabulation example
#add_custom_target(tabulation_example ALL
#  COMMAND ${Python3_EXECUTABLE} setup.py build_ext --inplace
#  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
#  DEPENDS ${SWIG_MODULE_hoppet_v1_REAL_NAME}
#)