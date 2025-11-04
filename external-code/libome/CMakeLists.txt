cmake_minimum_required(VERSION 3.20)
# Avoid warnings from FetchContent about timestamps
if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
  cmake_policy(SET CMP0135 NEW)
endif()

project(libome
  VERSION 1.0.2
  DESCRIPTION "Implementation of massive operator matrix elements (OMEs) of the QCD twist-2 operators in x-space"
)
set(CMAKE_CXX_STANDARD 17)
include(GNUInstallDirs)
option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

# Select "Release" build by default to enable optimisation
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Set the build type (Debug, Release, MinSizeRel, RelWithDebInfo)" FORCE)
endif()

# External dependencies
find_package(GSL REQUIRED)

# Floating point exceptions (FPEs)
# Note: enabling this changes nothing in the code, but if you want to rely on
# FPEs in your code, we have to make sure that the compiler does not make
# optimisations that would raise FPEs that the original code would not have
# raised. So far, this has turned out to be an issue with clang when compiling
# for Apple arm64 CPUs. See
# https://clang.llvm.org/docs/UsersManual.html#controlling-floating-point-behavior
option(ENABLE_FPES "Enable floating point exceptions")
if(ENABLE_FPES)
  add_compile_options("$<$<COMPILE_LANG_AND_ID:CXX,Clang>:-ffp-exception-behavior=maytrap>")
  add_compile_options("$<$<COMPILE_LANG_AND_ID:CXX,AppleClang>:-ffp-exception-behavior=maytrap>")
endif()

# Unit and integration test suite
option(BUILD_TESTS "Build the test suite" OFF)
if(BUILD_TESTS)
  include(FetchContent)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.17.0.zip
    )
  # For Windows: Prevent overriding the parent project's compiler/linker settings
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  set(INSTALL_GTEST OFF)
  set(BUILD_GMOCK OFF)
  FetchContent_MakeAvailable(googletest)
  
  enable_testing()
  include(GoogleTest)
  
  function(add_test_executable name)
    add_executable(${name} ${ARGN})
    target_link_libraries(${name} PRIVATE ome GTest::gtest_main)
    gtest_discover_tests(${name})
  endfunction()
  
  add_subdirectory(tests)
endif()

# Main code tree
add_subdirectory(src)

# Usage examples
option(BUILD_EXAMPLES "Build the example programs" OFF)
if(BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

# Provide pkg-config information for this library
cmake_path(APPEND exec_prefix_for_pc_file "\${prefix}" "${CMAKE_INSTALL_BINDIR}")
cmake_path(APPEND libdir_for_pc_file "\${prefix}" "${CMAKE_INSTALL_LIBDIR}")
cmake_path(APPEND includedir_for_pc_file "\${prefix}" "${CMAKE_INSTALL_INCLUDEDIR}")
configure_file(pkg-config.in ${CMAKE_PROJECT_NAME}.pc @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/${CMAKE_PROJECT_NAME}.pc
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
