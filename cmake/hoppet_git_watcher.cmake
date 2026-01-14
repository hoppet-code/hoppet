# This file derives from the original
#
# git_watcher.cmake
# https://raw.githubusercontent.com/andrew-hardin/cmake-git-version-tracking/master/git_watcher.cmake
#
# and then substantially refactored wtih ChatGPT 5.2 to adapt to modern CMake
# practices and to provide a cleaner public API, including better support
# for re-use across other projects.

# Originally released under the MIT License
# https://raw.githubusercontent.com/andrew-hardin/cmake-git-version-tracking/master/LICENSE
#
# Still released under the MIT License


cmake_minimum_required(VERSION 3.15)

# =============================================================================
# Public API
#
# git_watcher_add(
#   TARGET            <target-name>         # required, name of the custom target
#   PRE_CONFIGURE     <input-file>          # required, relative to current source
#   POST_CONFIGURE    <output-file>         # required, relative to current binary
#   [STATE_FILE       <file>]
#   [WORKING_DIR      <dir>]
#   [GIT_EXECUTABLE   <path>]
#   [FAIL_ON_ERROR    <ON|OFF>]
#   [VERBOSE          <ON|OFF>]             # outputs some status messages (default OFF)
#   [ENABLE           <ON|OFF>]             # whether to enable the target (default ON)
# )
# =============================================================================

# needed so that the function below can correctly identify the directory this file is
# in 
set(HOPPET_GW_MODULE_FILE "${CMAKE_CURRENT_LIST_FILE}")


function(hoppet_git_watcher_add)
  set(oneValueArgs
    TARGET
    PRE_CONFIGURE
    POST_CONFIGURE
    STATE_FILE
    WORKING_DIR
    GIT_EXECUTABLE
    FAIL_ON_ERROR
    VERBOSE
    ENABLE
  )
  cmake_parse_arguments(GW "" "${oneValueArgs}" "" ${ARGN})

  foreach(req TARGET PRE_CONFIGURE POST_CONFIGURE)
    if(NOT GW_${req})
      message(FATAL_ERROR "hoppet_git_watcher_add(): ${req} is required")
    endif()
  endforeach()

  get_filename_component(GW_PRE_CONFIGURE "${GW_PRE_CONFIGURE}"
    ABSOLUTE  BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
  )
  
  get_filename_component(GW_POST_CONFIGURE "${GW_POST_CONFIGURE}"
    ABSOLUTE  BASE_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  )

  if(NOT DEFINED GW_ENABLE)
    set(GW_ENABLE ON)
  endif()

  if (NOT GW_ENABLE)
    configure_file("${GW_PRE_CONFIGURE}" "${GW_POST_CONFIGURE}" COPYONLY)
    return()
  endif()

  if(TARGET ${GW_TARGET})
    # idempotent
    return()
  endif()

  if(NOT GW_WORKING_DIR)
    set(GW_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
  endif()

  if(NOT GW_STATE_FILE)
    set(GW_STATE_FILE "${CMAKE_BINARY_DIR}/${GW_TARGET}-git-state")
  endif()

  if(NOT DEFINED GW_FAIL_ON_ERROR)
    set(GW_FAIL_ON_ERROR ON)
  endif()

  if(NOT DEFINED GW_VERBOSE)
    set(GW_VERBOSE OFF)
  endif()

  if(NOT GW_GIT_EXECUTABLE)
    find_package(Git QUIET)
    set(GW_GIT_EXECUTABLE "${GIT_EXECUTABLE}")
  endif()

  if(NOT GW_GIT_EXECUTABLE)
    if(GW_FAIL_ON_ERROR)
      message(FATAL_ERROR "git executable not found")
    else()
      configure_file("${GW_PRE_CONFIGURE}" "${GW_POST_CONFIGURE}" COPYONLY)
      return()
    endif()
  endif()
  
  # this should work in older CMake versions
  set(_GW_SCRIPT "${HOPPET_GW_MODULE_FILE}")
  # this needs CMake 3.17
  #set(_GW_SCRIPT "${CMAKE_CURRENT_FUNCTION_LIST_FILE}")

  if (GW_VERBOSE)
    message(STATUS "git_watcher_add() for ${GW_TARGET}: GW_PRE_CONFIGURE='${GW_PRE_CONFIGURE}'")
    message(STATUS "git_watcher_add() for ${GW_TARGET}: GW_POST_CONFIGURE='${GW_POST_CONFIGURE}'")
  endif()

  add_custom_target(${GW_TARGET}
    ALL
    BYPRODUCTS
      "${GW_POST_CONFIGURE}"
      "${GW_STATE_FILE}"
    COMMENT "Checking git repository for changes (${GW_TARGET})"
    COMMAND
      ${CMAKE_COMMAND}
      -DGW_BUILD_TIME=TRUE
      -DGW_PRE_CONFIGURE=${GW_PRE_CONFIGURE}
      -DGW_POST_CONFIGURE=${GW_POST_CONFIGURE}
      -DGW_STATE_FILE=${GW_STATE_FILE}
      -DGW_WORKING_DIR=${GW_WORKING_DIR}
      -DGW_GIT_EXECUTABLE=${GW_GIT_EXECUTABLE}
      -DGW_FAIL_ON_ERROR=${GW_FAIL_ON_ERROR}
      -P "${_GW_SCRIPT}"
  )
endfunction()

# =============================================================================
# Build-time implementation (executed when cmake -P is invoked)
# =============================================================================

if(GW_BUILD_TIME)

  # Helper: run git command, return rc and output
  function(_gw_run_git out_var)
    execute_process(
      COMMAND "${GW_GIT_EXECUTABLE}" ${ARGN}
      WORKING_DIRECTORY "${GW_WORKING_DIR}"
      RESULT_VARIABLE _gw_rc
      OUTPUT_VARIABLE _gw_out
      ERROR_VARIABLE _gw_err
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(NOT _gw_rc EQUAL 0)
      set(GIT_RETRIEVED_STATE "false")
      if(GW_FAIL_ON_ERROR)
        # Show the stderr in a fatal error with the invoked args
        string(REPLACE ";" " " _gw_args "${ARGV}")
        message(FATAL_ERROR "${_gw_err} (${GW_GIT_EXECUTABLE} ${_gw_args})")
      endif()
    endif()
    set(${out_var} "${_gw_out}" PARENT_SCOPE)
  endfunction()

  # GetGitState: populate the GIT_* variables exactly like original
  function(_gw_get_git_state)
    # assume retrieved until an error occurs
    set(GIT_RETRIEVED_STATE "true" PARENT_SCOPE) # will be overwritten if error

    # default working dir already set via GW_WORKING_DIR
    # 1) status --porcelain -uno
    _gw_run_git(_gw_status_output status --porcelain -uno)
    if(NOT DEFINED _gw_rc) # _gw_rc only used for fatal above; check output variable instead
      # proceed
    endif()

    if(NOT "${_gw_status_output}" STREQUAL "")
      # replace CR/LF with semicolon as in original
      string(REGEX REPLACE "[\r\n]+" ";" _gw_status_proc "${_gw_status_output}")
      set(GIT_STATUS_UNO "${_gw_status_proc}" PARENT_SCOPE)
      set(GIT_IS_DIRTY "true" PARENT_SCOPE)

      # create GIT_F90_STATUS_UNO similar to original:
      # the original produced something like: replace ";" with ;"//&\n     "
      # Reproduce same replacement string (including quoting)
      string(REPLACE ";" ";\"//&\n     \"" _gw_f90 "${_gw_status_proc}")
      set(GIT_F90_STATUS_UNO "${_gw_f90}" PARENT_SCOPE)
    else()
      set(GIT_STATUS_UNO "" PARENT_SCOPE)
      set(GIT_F90_STATUS_UNO "" PARENT_SCOPE)
      set(GIT_IS_DIRTY "false" PARENT_SCOPE)
    endif()

    # object = HEAD (original)
    set(_gw_object "HEAD")

    # HEAD SHA1
    _gw_run_git(_gw_out show -s --format=%H ${_gw_object})
    if(DEFINED _gw_out)
      set(GIT_HEAD_SHA1 "${_gw_out}" PARENT_SCOPE)
    endif()

    # author name
    _gw_run_git(_gw_out show -s --format=%an ${_gw_object})
    if(DEFINED _gw_out)
      set(GIT_AUTHOR_NAME "${_gw_out}" PARENT_SCOPE)
    endif()

    # author email
    _gw_run_git(_gw_out show -s --format=%ae ${_gw_object})
    if(DEFINED _gw_out)
      set(GIT_AUTHOR_EMAIL "${_gw_out}" PARENT_SCOPE)
    endif()

    # commit date ISO8601
    _gw_run_git(_gw_out show -s --format=%ci ${_gw_object})
    if(DEFINED _gw_out)
      set(GIT_COMMIT_DATE_ISO8601 "${_gw_out}" PARENT_SCOPE)
    endif()

    # commit subject (escape quotes)
    _gw_run_git(_gw_out show -s --format=%s ${_gw_object})
    if(DEFINED _gw_out)
      string(REPLACE "\"" "\\\\\"" _gw_escaped "${_gw_out}")
      set(GIT_COMMIT_SUBJECT "${_gw_escaped}" PARENT_SCOPE)
    endif()

    # commit body (escape quotes and escape line breaks similarly to original)
    _gw_run_git(_gw_out show -s --format=%b ${_gw_object})
    if(DEFINED _gw_out AND NOT _gw_out STREQUAL "")
      # Escape double quotes
      string(REPLACE "\"" "\\\\\"" _gw_body_escaped "${_gw_out}")

      # Escape CRLF first; then LF if needed, replicating original behavior
      string(REPLACE "\r\n" "\\r\\n\\\r\n" _gw_body_safe "${_gw_body_escaped}")
      if(_gw_body_safe STREQUAL _gw_body_escaped)
        # no CRLF changes, try LF
        string(REPLACE "\n" "\\n\\\n" _gw_body_safe "${_gw_body_escaped}")
      endif()
      set(GIT_COMMIT_BODY "\"${_gw_body_safe}\"" PARENT_SCOPE)
    else()
      set(GIT_COMMIT_BODY "\"\"" PARENT_SCOPE)
    endif()

    # git describe --always
    _gw_run_git(_gw_out describe --always ${_gw_object})
    if(NOT DEFINED _gw_out OR _gw_out STREQUAL "")
      set(GIT_DESCRIBE "unknown" PARENT_SCOPE)
    else()
      set(GIT_DESCRIBE "${_gw_out}" PARENT_SCOPE)
    endif()

    # mark retrieved true unless an earlier RunGitCommand fataled
    if(NOT DEFINED GIT_RETRIEVED_STATE)
      set(GIT_RETRIEVED_STATE "true" PARENT_SCOPE)
    endif()
  endfunction()

  # HashGitState: compute hash identical to original algorithm
  function(_gw_hash_git_state out_var)
    # list of names exactly as original order
    set(_gw_state_names
      GIT_RETRIEVED_STATE
      GIT_HEAD_SHA1
      GIT_IS_DIRTY
      GIT_STATUS_UNO
      GIT_F90_STATUS_UNO
      GIT_AUTHOR_NAME
      GIT_AUTHOR_EMAIL
      GIT_COMMIT_DATE_ISO8601
      GIT_COMMIT_SUBJECT
      GIT_COMMIT_BODY
      GIT_DESCRIBE
    )

    set(_gw_acc "")
    foreach(_gw_n ${_gw_state_names})
      # avoid undefined-variable expansion producing errors: replace with empty string if not defined
      if(DEFINED ${_gw_n})
        set(_gw_val "${${_gw_n}}")
      else()
        set(_gw_val "")
      endif()
      string(SHA256 _gw_acc "${_gw_acc}${_gw_val}")
      # move computed to the accumulator variable
      set(_gw_acc "${_gw_acc}")
    endforeach()
    set(${out_var} "${_gw_acc}" PARENT_SCOPE)
  endfunction()

  # CheckGit: check if state changed and update state file + configure
  function(_gw_check_and_update)
    # populate GIT_* variables
    _gw_get_git_state()

    # compute state hash
    _gw_hash_git_state(_gw_state)

    # include PRE_CONFIGURE file hash into state (exact parity with original)
    file(SHA256 "${GW_PRE_CONFIGURE}" _gw_pre_hash)
    string(SHA256 _gw_state "${_gw_pre_hash}${_gw_state}")

    # compare against existing file
    if(EXISTS "${GW_STATE_FILE}")
      file(READ "${GW_STATE_FILE}" _gw_old)
      if(_gw_old STREQUAL "${_gw_state}")
        # no change
        set(_gw_changed "false" PARENT_SCOPE)
        return()
      endif()
    endif()

    # state changed -> write file and configure_file
    file(WRITE "${GW_STATE_FILE}" "${_gw_state}")
    set(_gw_changed "true" PARENT_SCOPE)

    # Before configure_file, expose variables in this script scope exactly as original expects:
    # The original set variables (not env) before calling configure_file; do same:
    # GIT_RETRIEVED_STATE, GIT_HEAD_SHA1, etc are already set in this scope by _gw_get_git_state()

    configure_file("${GW_PRE_CONFIGURE}" "${GW_POST_CONFIGURE}" @ONLY)
  endfunction()

  # Main build-time flow (parity)
  # if the post-configure file does not exist, we must regenerate as original did
  if(NOT EXISTS "${GW_POST_CONFIGURE}")
    # run check and update; it will create POST and state file
    _gw_check_and_update()
  else()
    # run check; only updates if changed
    _gw_check_and_update()
  endif()

endif() # GW_BUILD_TIME
