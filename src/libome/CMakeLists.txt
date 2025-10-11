set(LIBOME_DIR ${PROJECT_SOURCE_DIR}/submodules/libome-fork/src)
set(CXXINCLUDE ${CXXINCLUDE} ${LIBOME_DIR} PARENT_SCOPE)

file(GLOB LIBOME_SOURCES_CXX ${LIBOME_DIR}/ome/A*.cpp)
set(LIBOME_SOURCES_CXX ${LIBOME_SOURCES_CXX} ${PROJECT_SOURCE_DIR}/src/libome/hoppet_libome_interface.cc)

# make sure we modify the parent-scope variable
set(SOURCESCXX ${SOURCESCXX} ${LIBOME_SOURCES_CXX} PARENT_SCOPE)
set(SOURCESf90 ${SOURCESf90} src/libome/hoppet_libome_interfaces.f90 PARENT_SCOPE)