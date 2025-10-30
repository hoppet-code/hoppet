
# Include(FetchContent)
# FetchContent_Declare(
#   Catch2
#   GIT_REPOSITORY https://github.com/catchorg/Catch2.git
#   GIT_TAG        v3.8.1 # or a later release
# )
# FetchContent_MakeAvailable(Catch2)

## Attempt to find Catch2; if not found, skip building unit tests.
## Install with `brew install catch2` or `apt-get install libcatch2-dev`.
find_package(Catch2 3 QUIET)

if(Catch2_FOUND)
	set(UNIT_TESTS_SOURCES unit-tests.cc)
	add_executable(hoppet-unit-tests-cxx ${UNIT_TESTS_SOURCES})
	target_include_directories(hoppet-unit-tests-cxx PUBLIC ${PROJECT_SOURCE_DIR}/src)
	target_link_libraries(hoppet-unit-tests-cxx PRIVATE Catch2::Catch2WithMain)
	target_link_libraries(hoppet-unit-tests-cxx PRIVATE hoppet_static)
	add_test(NAME hoppet-unit-tests-cxx COMMAND hoppet-unit-tests-cxx)
else()
	message(STATUS "Catch2 not found - skipping C++ unit tests (unit-tests-cxx)")
endif()