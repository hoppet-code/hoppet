cmake_minimum_required(VERSION 3.16...3.31)

set(UseSWIG_MODULE_VERSION 2)
include(UseSWIG)

if (CMAKE_VERSION VERSION_LESS 3.18)
  set(DEV_MODULE Development)
else()
  set(DEV_MODULE Development.Module)
endif()

find_package(Python 3.8 COMPONENTS Interpreter ${DEV_MODULE} REQUIRED)

set(NUMPY_SWIG_DIR "${PROJECT_SOURCE_DIR}/pyinterface/inc")
set_property(SOURCE hoppet.i PROPERTY SWIG_FLAGS "-I${NUMPY_SWIG_DIR}")

set_property(SOURCE hoppet.i PROPERTY CPLUSPLUS ON)
# Remove or set to your actual module name, not __init__
set_property(SOURCE hoppet.i PROPERTY SWIG_MODULE_NAME __init__)

swig_add_library(hoppet_swig LANGUAGE python SOURCES hoppet.i OUTPUT_DIR hoppet_swig)
set_property(TARGET hoppet_swig PROPERTY SWIG_USE_TARGET_INCLUDE_DIRECTORIES TRUE)
target_include_directories(hoppet_swig PRIVATE ${NUMPY_SWIG_DIR} ${PROJECT_SOURCE_DIR}/src)

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set_target_properties(hoppet_swig PROPERTIES INSTALL_RPATH "@loader_path/../lib;@loader_path/../../../../lib")
else()
  set_target_properties(hoppet_swig PROPERTIES INSTALL_RPATH "$ORIGIN/../lib:$ORIGIN/../../../../lib:$ORIGIN/../../../../lib64")
endif()

set_property(TARGET hoppet_swig PROPERTY SUFFIX ".${Python_SOABI}${CMAKE_SHARED_MODULE_SUFFIX}")

#target_link_libraries(hoppet_swig PUBLIC hoppet)
target_link_libraries(hoppet_swig PUBLIC hoppet_static)
target_link_libraries(hoppet_swig PUBLIC Python::Module)

if ("${SWIG_VERSION}" VERSION_LESS "4.2.0")
  target_compile_definitions(hoppet_swig PUBLIC -DSWIG)
  target_compile_definitions(hoppet_swig PUBLIC -DSWIG_TYPE_TABLE=hoppet_swig)
endif()

target_include_directories(hoppet_swig 
                           PUBLIC 
                            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                            $<INSTALL_INTERFACE:include>
                           )
add_library(hoppet::hoppet_swig ALIAS hoppet_swig)

if (HOPPET_USE_PYTHON_SITEARCH)
  set(HOPPET_PYTHON_INSTALL_DIR "${Python_SITEARCH}")
elseif(NOT "${HOPPET_CUSTOM_PYTHON_INSTALL}" STREQUAL "")
  set(HOPPET_PYTHON_INSTALL_DIR "${HOPPET_CUSTOM_PYTHON_INSTALL}")
else()
  set(HOPPET_PYTHON_INSTALL_DIR "lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages")
endif()

set(pyexecdir "${HOPPET_PYTHON_INSTALL_DIR}" PARENT_SCOPE) # PARENT_SCOPE tells CMake to pass this back to the parent CMakeLists.txt
if(HOPPET_BUILD_PYINTERFACE)
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/hoppet_swig/__init__.py" DESTINATION "${HOPPET_PYTHON_INSTALL_DIR}/${HOPPET_PYTHON_PACKAGE_NAME}")
  install(FILES "${PROJECT_SOURCE_DIR}/pyinterface/lhapdf.py" DESTINATION "${HOPPET_PYTHON_INSTALL_DIR}/${HOPPET_PYTHON_PACKAGE_NAME}")
  install(TARGETS hoppet_swig EXPORT HoppetTargets DESTINATION "${HOPPET_PYTHON_INSTALL_DIR}/${HOPPET_PYTHON_PACKAGE_NAME}")
elseif(HOPPET_BUILD_PYWHEELS)
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/hoppet_swig/__init__.py" DESTINATION hoppet)
  install(FILES "${PROJECT_SOURCE_DIR}/pyinterface/lhapdf.py" DESTINATION hoppet)
  install(TARGETS hoppet_swig EXPORT HoppetTargets DESTINATION hoppet)
endif()
  #install(TARGETS hoppet DESTINATION hoppet)
set(HOPPET_PYTHON_INSTALL_DIR ${HOPPET_PYTHON_INSTALL_DIR} PARENT_SCOPE) 

