add_executable(hoppet-unit-tests unit_tests.f90 hoppet-unit-tests.f90)
target_link_libraries(hoppet-unit-tests hoppet)
add_test(NAME hoppet-unit-tests COMMAND hoppet-unit-tests)