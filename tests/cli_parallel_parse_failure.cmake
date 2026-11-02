# A malformed second gene tree, parsed with two workers, must end the run
# with the handled exit status 1 and report the parser error.

if(NOT DEFINED QRAFT_CLI OR NOT DEFINED QRAFT_CLI_TEST_DIR)
  message(FATAL_ERROR "QRAFT_CLI and QRAFT_CLI_TEST_DIR are required")
endif()

set(input "${QRAFT_CLI_TEST_DIR}/parallel_parse_failure.genes.nwk")
set(output "${QRAFT_CLI_TEST_DIR}/parallel_parse_failure.species.nwk")
file(WRITE "${input}" "((A,B),(C,D));\n((A,B),(C,D);\n")
file(REMOVE "${output}")

execute_process(
  COMMAND "${QRAFT_CLI}"
          --genes "${input}"
          --out "${output}"
          --threads 2
          --no-progress
  RESULT_VARIABLE result
  OUTPUT_VARIABLE stdout
  ERROR_VARIABLE stderr)

file(REMOVE "${input}" "${output}")

if(NOT result STREQUAL "1")
  message(FATAL_ERROR
    "parallel parse failure returned '${result}', expected handled exit 1; "
    "stdout='${stdout}' stderr='${stderr}'")
endif()
if(NOT stderr MATCHES "expected")
  message(FATAL_ERROR
    "parallel parse failure did not report the parser error; stderr='${stderr}'")
endif()
