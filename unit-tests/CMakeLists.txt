add_executable(hoppet-unit-tests unit_tests.f90 hoppet-unit-tests.F90)
target_link_libraries(hoppet-unit-tests hoppet)
add_test(NAME hoppet-unit-tests COMMAND hoppet-unit-tests)

add_executable(check-thread-safety check-thread-safety.cc)
target_include_directories(check-thread-safety PUBLIC ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(check-thread-safety hoppet)
