
set(UNIT_TESTS_SOURCES unit-tests.cc)
add_executable(hoppet-unit-tests-cxx ${UNIT_TESTS_SOURCES})
target_include_directories(hoppet-unit-tests-cxx PUBLIC ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(hoppet-unit-tests-cxx hoppet_static)
add_test(NAME hoppet-unit-tests-cxx COMMAND hoppet-unit-tests-cxx)