add_executable(sample-omes sample-omes.cpp)
target_link_libraries(sample-omes PUBLIC ome)

add_executable(access-coefficients access-coefficients.cpp)
target_link_libraries(access-coefficients PUBLIC ome)

add_executable(calculate-moments calculate-moments.cpp)
target_link_libraries(calculate-moments PUBLIC ome)

add_executable(calculate-convolution calculate-convolution.cpp)
target_link_libraries(calculate-convolution PUBLIC ome)

add_executable(use-c-interface use-c-interface.c)
target_link_libraries(use-c-interface PUBLIC ome m)

option(BUILD_FORTRAN_EXAMPLE "Build the example how to use the libarary from Fortran" OFF)
if(BUILD_FORTRAN_EXAMPLE)
  enable_language(Fortran)
  add_executable(use-c-interface-from-fortran use-c-interface-from-fortran.f90)
  target_link_libraries(use-c-interface-from-fortran PUBLIC ome)
endif()
