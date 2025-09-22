set(HELPERS unit_tests.f90 test_alphas.f90 test_interpolation.F90 
            ../benchmarking/io_utils.f90)
add_executable(hoppet-unit-tests hoppet-unit-tests.F90 ${HELPERS})
target_link_libraries(hoppet-unit-tests hoppet)
add_test(NAME hoppet-unit-tests COMMAND hoppet-unit-tests)