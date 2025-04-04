cmake_minimum_required(VERSION 3.15)
project(HoppetPython)

#------------------------------------------------------------------------
# basic requirements
SET(CMAKE_CXX_STANDARD 14)
SET(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "-O2 -g -Wall")


# specify that we're using fortran
enable_language(Fortran)
set(CMAKE_Fortran_FLAGS "-O3")


SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

#------------------------------------------------------------------------
# find Swig
#
# The following lines get rid of warnings
# on some systems (Gavin's and Gregory's swig 4)
if(POLICY CMP0078)
  cmake_policy(SET CMP0078 NEW)
endif()
if(POLICY CMP0086)
  cmake_policy(SET CMP0086 NEW)
endif()

FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

#------------------------------------------------------------------------
# Find Python3
#
# the more modern way (as of 2021) of finding Python3 components
# (for Gavin this worked out of the box, whereas the previous
# solution needed manual specification of non-standard lib and
# include paths)
find_package (Python3 COMPONENTS Interpreter Development REQUIRED)
message("Python3 lib: " ${Python3_LIBRARIES})
message("Python3 include: " ${Python3_INCLUDE_DIRS})
message("Python3 site lib: " ${Python3_SITELIB})
INCLUDE_DIRECTORIES(${Python3_INCLUDE_DIRS})

# Note: relative paths are all given wrt to the sourcedir
set(PARENT_DIR .. CACHE PATH "Path to the main (hoppet) directory")
message ("-- Looking for dependences in parent directory ${PARENT_DIR}")


#----------------------------------------------------------------------
# find the Hoppet directories
#
#   HOPPET_SOURCE_DIR  where most includes will be found [${PARENT_DIR}/src]
#   HOPPET_BUILD_DIR   where the library is found [${PARENT_DIR}/build/lib]

# find the source dir
find_path( HOPPET_SOURCE_DIR 
  NAMES hoppet_v1.h
  HINTS ${PARENT_DIR}/src
)
if (HOPPET_SOURCE_DIR)
  message ("-- Hoppet sources (includes) in ${HOPPET_SOURCE_DIR}")
else()
  message(FATAL_ERROR "Hoppet sources not found")
endif()

# find the include dirs for the config.hh file
find_path( HOPPET_BUILD_DIR
  NAMES bin
  HINTS ${PARENT_DIR}/build
)
if (HOPPET_BUILD_DIR)
  message ("-- Hoppet build in ${HOPPET_BUILD_DIR}")
else()
  message(FATAL_ERROR "Hoppet build not found")
endif()

# find the Hoppet library
find_library( HOPPET_LIBRARIES 
  NAMES hoppet_v1
  HINTS ${PARENT_BUILD_DIR}
)
if (HOPPET_LIBRARIES)
  message("-- Hoppet library ${HOPPET_LIBRARIES}")
else()
  message(FATAL_ERROR "Hoppet library not found")
endif()

#------------------------------------------------------------------------
INCLUDE_DIRECTORIES(${HOPPET_SOURCE_DIR} ${HOPPET_BUILD_DIR})
SET_SOURCE_FILES_PROPERTIES(hoppet_v1.i PROPERTIES CPLUSPLUS ON)
SWIG_ADD_LIBRARY(hoppet_v1 LANGUAGE python SOURCES hoppet_v1.i)
SWIG_LINK_LIBRARIES(hoppet_v1 ${Python3_LIBRARIES} ${HOPPET_LIBRARIES})

message("-- SWIG module name ${SWIG_MODULE_hoppet_v1_REAL_NAME}")

install(TARGETS ${SWIG_MODULE_hoppet_v1_REAL_NAME}
  LIBRARY DESTINATION ${Python3_SITEARCH}
  )

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/hoppet_v1.py 
          DESTINATION ${Python3_SITELIB})
