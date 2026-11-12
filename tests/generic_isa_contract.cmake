# A generic build must carry no -march=native, -mcpu=native or -mtune=native
# flag; the probes first check that the matcher recognizes each spelling. A
# GCC or clang build must compile the library with -ffp-contract=off. A
# jemalloc-linked library must not reference glibc's malloc_trim or mallopt.

if(NOT DEFINED QRAFT_FLAG_MANIFEST OR
   NOT EXISTS "${QRAFT_FLAG_MANIFEST}")
  message(FATAL_ERROR "missing generated compile-flag manifest")
endif()

file(READ "${QRAFT_FLAG_MANIFEST}" QRAFT_EFFECTIVE_FLAGS)
set(QRAFT_NATIVE_CPU_FLAG_PATTERN
  "(^|[^A-Za-z0-9_])-(march|mcpu|mtune)=native($|[^A-Za-z0-9_])")
foreach(QRAFT_MATCHER_PROBE IN ITEMS
    "global=-march=native"
    "library=SHELL:-mcpu=native "
    "link=LINKER:-mtune=native,")
  if(NOT QRAFT_MATCHER_PROBE MATCHES "${QRAFT_NATIVE_CPU_FLAG_PATTERN}")
    message(FATAL_ERROR "native CPU flag matcher missed one of its probes")
  endif()
endforeach()

if(QRAFT_EFFECTIVE_FLAGS MATCHES "${QRAFT_NATIVE_CPU_FLAG_PATTERN}")
  message(FATAL_ERROR
    "generic build contains machine-local CPU tuning:\n"
    "${QRAFT_EFFECTIVE_FLAGS}")
endif()

# Builds with GCC or clang compile without floating-point contraction, so a
# fused multiply-add cannot change a score on aarch64 or Apple silicon.
if(QRAFT_EXPECT_FP_CONTRACT_OFF)
  string(REGEX MATCH "(^|\n)library=[^\n]*" QRAFT_LIBRARY_FLAGS
    "${QRAFT_EFFECTIVE_FLAGS}")
  if(NOT QRAFT_LIBRARY_FLAGS MATCHES "-ffp-contract=off")
    message(FATAL_ERROR
      "the library is not compiled with -ffp-contract=off:\n"
      "${QRAFT_EFFECTIVE_FLAGS}")
  endif()
endif()

if(QRAFT_EFFECTIVE_FLAGS MATCHES "QRAFT_LINKED_JEMALLOC(=1)?")
  if(NOT DEFINED QRAFT_LIBRARY OR NOT EXISTS "${QRAFT_LIBRARY}")
    message(FATAL_ERROR
      "jemalloc-linked build is missing the project library archive")
  endif()
  if(NOT DEFINED QRAFT_NM OR NOT EXISTS "${QRAFT_NM}")
    message(FATAL_ERROR
      "jemalloc-linked build requires nm for its allocator-symbol contract")
  endif()
  execute_process(
    COMMAND "${QRAFT_NM}" -u "${QRAFT_LIBRARY}"
    RESULT_VARIABLE QRAFT_NM_RESULT
    OUTPUT_VARIABLE QRAFT_UNDEFINED_SYMBOLS
    ERROR_VARIABLE QRAFT_NM_ERROR)
  if(NOT QRAFT_NM_RESULT EQUAL 0)
    message(FATAL_ERROR
      "nm could not inspect the jemalloc-linked project archive:\n"
      "${QRAFT_NM_ERROR}")
  endif()
  if(QRAFT_UNDEFINED_SYMBOLS MATCHES
     "(^|[\r\n])[^\r\n]*(malloc_trim|mallopt)($|[\r\n])")
    message(FATAL_ERROR
      "jemalloc-linked code still references glibc allocator controls:\n"
      "${QRAFT_UNDEFINED_SYMBOLS}")
  endif()
endif()
