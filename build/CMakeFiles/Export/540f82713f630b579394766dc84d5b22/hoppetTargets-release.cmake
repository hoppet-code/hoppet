#----------------------------------------------------------------
# Generated CMake target import file for configuration "RELEASE".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "hoppet::hoppet" for configuration "RELEASE"
set_property(TARGET hoppet::hoppet APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hoppet::hoppet PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libhoppet.so.0.0.0"
  IMPORTED_SONAME_RELEASE "libhoppet.so.0.0.0"
  )

list(APPEND _cmake_import_check_targets hoppet::hoppet )
list(APPEND _cmake_import_check_files_for_hoppet::hoppet "${_IMPORT_PREFIX}/lib/libhoppet.so.0.0.0" )

# Import target "hoppet::hoppet_static" for configuration "RELEASE"
set_property(TARGET hoppet::hoppet_static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hoppet::hoppet_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "Fortran"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libhoppet.a"
  )

list(APPEND _cmake_import_check_targets hoppet::hoppet_static )
list(APPEND _cmake_import_check_files_for_hoppet::hoppet_static "${_IMPORT_PREFIX}/lib/libhoppet.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
