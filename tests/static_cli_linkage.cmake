# The static CLI must have no PT_INTERP segment and no DT_NEEDED entry,
# must contain jemalloc's mallctl when jemalloc is expected, and must run.

set(ENV{LC_ALL} "C")

if(NOT DEFINED QRAFT_CLI OR NOT EXISTS "${QRAFT_CLI}")
  message(FATAL_ERROR "missing qraft executable")
endif()
if(NOT DEFINED QRAFT_READELF OR NOT EXISTS "${QRAFT_READELF}")
  message(FATAL_ERROR "readelf is required for the static-linkage contract")
endif()

execute_process(
  COMMAND "${QRAFT_READELF}" -lW "${QRAFT_CLI}"
  RESULT_VARIABLE QRAFT_PROGRAM_HEADERS_RESULT
  OUTPUT_VARIABLE QRAFT_PROGRAM_HEADERS
  ERROR_VARIABLE QRAFT_PROGRAM_HEADERS_ERROR)
if(NOT QRAFT_PROGRAM_HEADERS_RESULT EQUAL 0)
  message(FATAL_ERROR
    "readelf could not inspect program headers: "
    "${QRAFT_PROGRAM_HEADERS_ERROR}")
endif()
if(QRAFT_PROGRAM_HEADERS MATCHES
   "(^|[\r\n])[^\r\n]*INTERP[^\r\n]*($|[\r\n])")
  message(FATAL_ERROR "static executable contains a PT_INTERP segment")
endif()

execute_process(
  COMMAND "${QRAFT_READELF}" -dW "${QRAFT_CLI}"
  RESULT_VARIABLE QRAFT_DYNAMIC_RESULT
  OUTPUT_VARIABLE QRAFT_DYNAMIC_SECTION
  ERROR_VARIABLE QRAFT_DYNAMIC_ERROR)
if(NOT QRAFT_DYNAMIC_RESULT EQUAL 0)
  message(FATAL_ERROR
    "readelf could not inspect the dynamic section: ${QRAFT_DYNAMIC_ERROR}")
endif()
if(QRAFT_DYNAMIC_SECTION MATCHES "\\(NEEDED\\)")
  message(FATAL_ERROR "static executable contains a DT_NEEDED dependency")
endif()

if(QRAFT_EXPECT_JEMALLOC)
  if(NOT DEFINED QRAFT_NM OR NOT EXISTS "${QRAFT_NM}")
    message(FATAL_ERROR
      "nm is required to prove jemalloc incorporation in the static CLI")
  endif()
  execute_process(
    COMMAND "${QRAFT_NM}" -g "${QRAFT_CLI}"
    RESULT_VARIABLE QRAFT_NM_RESULT
    OUTPUT_VARIABLE QRAFT_GLOBAL_SYMBOLS
    ERROR_VARIABLE QRAFT_NM_ERROR)
  if(NOT QRAFT_NM_RESULT EQUAL 0)
    message(FATAL_ERROR
      "nm could not inspect the static executable: ${QRAFT_NM_ERROR}")
  endif()
  if(NOT QRAFT_GLOBAL_SYMBOLS MATCHES
     "(^|[\r\n])[0-9A-Fa-f]+[ \t]+[A-Za-z][ \t]+mallctl($|[\r\n])")
    message(FATAL_ERROR
      "static executable does not contain jemalloc's mallctl symbol")
  endif()
endif()

if(NOT QRAFT_CROSSCOMPILING)
  execute_process(
    COMMAND "${QRAFT_CLI}" --help
    RESULT_VARIABLE QRAFT_HELP_RESULT
    OUTPUT_QUIET
    ERROR_QUIET)
  if(NOT QRAFT_HELP_RESULT EQUAL 0)
    message(FATAL_ERROR "static executable did not run successfully")
  endif()
endif()
