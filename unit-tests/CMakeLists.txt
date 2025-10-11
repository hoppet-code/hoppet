set(HELPERS unit_tests.f90 
    test_alphas.F90 
    test_interpolation.F90 
    test_convolution.F90
    test_convolution_C.cc
    ../benchmarking/io_utils.f90)
add_executable(hoppet-unit-tests hoppet-unit-tests.F90 ${HELPERS})
target_include_directories(hoppet-unit-tests PUBLIC ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(hoppet-unit-tests hoppet_static)
add_test(NAME hoppet-unit-tests COMMAND hoppet-unit-tests)

# Ensure the Fortran linker is used when this target contains C++ sources
# (so the Fortran runtime/main are properly linked). Without this CMake may
# choose the C++ linker and the Fortran main can be missing at link time.
set_property(TARGET hoppet-unit-tests PROPERTY LINKER_LANGUAGE Fortran)

add_executable(check-thread-safety check-thread-safety.cc)
target_include_directories(check-thread-safety PUBLIC ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(check-thread-safety hoppet_static)

# as of 2025-09-23, thread safety checks for ifx only work if
# it has the -recursive option (set in the make CMakeLists.txt)
add_test(NAME check-thread-safety COMMAND check-thread-safety)
