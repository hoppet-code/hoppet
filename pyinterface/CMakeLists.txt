cmake_minimum_required(VERSION 3.16...3.31)

if (CMAKE_VERSION VERSION_LESS 3.18)
  set(DEV_MODULE Development)
else()
  set(DEV_MODULE Development.Module)
endif()

find_package(Python 3.8 COMPONENTS Interpreter ${DEV_MODULE} REQUIRED)

# Find pybind11
find_package(pybind11 REQUIRED)

# Create the pybind11 module
pybind11_add_module(hoppet_pybind11 hoppet_pybind11.cpp)

# Set the output directory to create hoppet package
set_target_properties(hoppet_pybind11 PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/hoppet"
)

# Copy the __init__.py file
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/__init__.py ${CMAKE_CURRENT_BINARY_DIR}/hoppet/__init__.py COPYONLY)

# Include directories
target_include_directories(hoppet_pybind11 PRIVATE ${PROJECT_SOURCE_DIR}/src)

# NB this is specifically for building scikit-hep/hoppet wheels
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set_target_properties(hoppet_pybind11 PROPERTIES INSTALL_RPATH "@loader_path/../lib;@loader_path/../../../../lib")
else()
  set_target_properties(hoppet_pybind11 PROPERTIES INSTALL_RPATH "$ORIGIN/../lib:$ORIGIN/../../../../lib:$ORIGIN/../../../../lib64")
endif()

target_link_libraries(hoppet_pybind11 PUBLIC hoppet_static)

target_include_directories(hoppet_pybind11 
                           PUBLIC 
                            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                            $<INSTALL_INTERFACE:include>
                           )
add_library(hoppet::hoppet_pybind11 ALIAS hoppet_pybind11)

# Really should go to ${Python_SITEARCH} or into the wheel dir of scikit-build-core
if (HOPPET_USE_PYTHON_SITEARCH)
  set(HOPPET_PYTHON_INSTALL_DIR "${Python_SITEARCH}")
elseif(NOT "${HOPPET_CUSTOM_PYTHON_INSTALL}" STREQUAL "")
  set(HOPPET_PYTHON_INSTALL_DIR "${HOPPET_CUSTOM_PYTHON_INSTALL}")
else()
  set(HOPPET_PYTHON_INSTALL_DIR "lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages")
  #set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  #set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  #message(STATUS "CMAKE_INSTALL_RPATH: ${CMAKE_INSTALL_RPATH}")
endif()
set(pyexecdir "${HOPPET_PYTHON_INSTALL_DIR}" PARENT_SCOPE) # PARENT_SCOPE tells CMake to pass this back to the parent CMakeLists.txt
if(HOPPET_BUILD_PYINTERFACE)
  install(TARGETS hoppet_pybind11 EXPORT HoppetTargets DESTINATION "${HOPPET_PYTHON_INSTALL_DIR}/${HOPPET_PYTHON_PACKAGE_NAME}")
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/hoppet/__init__.py" DESTINATION "${HOPPET_PYTHON_INSTALL_DIR}/${HOPPET_PYTHON_PACKAGE_NAME}")
elseif(HOPPET_BUILD_PYWHEELS)
  install(TARGETS hoppet_pybind11 EXPORT HoppetTargets DESTINATION hoppet)
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/hoppet/__init__.py" DESTINATION hoppet)
endif()
  #install(TARGETS hoppet DESTINATION hoppet)
set(HOPPET_PYTHON_INSTALL_DIR ${HOPPET_PYTHON_INSTALL_DIR} PARENT_SCOPE) 

