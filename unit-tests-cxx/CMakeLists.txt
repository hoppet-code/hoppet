
# Include(FetchContent)
# FetchContent_Declare(
#   Catch2
#   GIT_REPOSITORY https://github.com/catchorg/Catch2.git
#   GIT_TAG        v3.8.1 # or a later release
# )
# FetchContent_MakeAvailable(Catch2)

## Attempt to find Catch2; if not found, fetch it via FetchContent.
## Install via package manager (brew/apt) to avoid fetching at configure time.
find_package(Catch2 3 QUIET)

if(NOT Catch2_FOUND)
	include(FetchContent)
	FetchContent_Declare(
		Catch2
		GIT_REPOSITORY https://github.com/catchorg/Catch2.git
		GIT_TAG        v3.11.0 # or a later release
	)
	# This will download and add Catch2 targets (Catch2::Catch2WithMain)
	FetchContent_MakeAvailable(Catch2)
endif()

if(Catch2_FOUND OR TARGET Catch2::Catch2WithMain)
	set(UNIT_TESTS_SOURCES unit-tests.cc)
	add_executable(hoppet-unit-tests-cxx ${UNIT_TESTS_SOURCES})
	target_include_directories(hoppet-unit-tests-cxx PUBLIC ${PROJECT_SOURCE_DIR}/src)
	target_link_libraries(hoppet-unit-tests-cxx PRIVATE Catch2::Catch2WithMain)
	target_link_libraries(hoppet-unit-tests-cxx PRIVATE hoppet_static)
	add_test(NAME hoppet-unit-tests-cxx COMMAND hoppet-unit-tests-cxx)
else()
	message(STATUS "Catch2 not found and could not be fetched - skipping C++ unit tests (unit-tests-cxx)")
endif()