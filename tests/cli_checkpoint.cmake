# Checks the decision journal (--checkpoint, --resume): a journaled run gives
# the same tree as a plain run, and a run killed right after, or in the middle
# of, a journal write and then resumed gives that tree too, at another thread
# count. The kills use the test hooks QRAFT_CHECKPOINT_TEST_KILL_AFTER and
# QRAFT_CHECKPOINT_TEST_KILL_MID, which end the process like kill -9. A due
# journal write is not held back by a long step.

if(NOT DEFINED QRAFT_CLI OR NOT DEFINED QRAFT_SOURCE_DIR OR
   NOT DEFINED QRAFT_CLI_TEST_DIR)
  message(FATAL_ERROR "QRAFT_CLI, QRAFT_SOURCE_DIR, and QRAFT_CLI_TEST_DIR are required")
endif()

set(input "${QRAFT_SOURCE_DIR}/tests/trees/astral2_200taxa_50genes.tre")
set(output_dir "${QRAFT_CLI_TEST_DIR}/checkpoint")
file(REMOVE_RECURSE "${output_dir}")
file(MAKE_DIRECTORY "${output_dir}")

# Runs the CLI with the extra arguments; `env` is a list of NAME=VALUE pairs.
function(run_qraft name expect_success env)
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E env ${env}
            "${QRAFT_CLI}" --genes "${input}" --no-progress --quiet ${ARGN}
    RESULT_VARIABLE result
    OUTPUT_VARIABLE stdout
    ERROR_VARIABLE stderr)
  if(expect_success AND NOT result EQUAL 0)
    message(FATAL_ERROR
      "${name} failed with ${result}\nstdout:\n${stdout}\nstderr:\n${stderr}")
  endif()
  if(NOT expect_success AND result EQUAL 0)
    message(FATAL_ERROR "${name} was expected to be killed but finished")
  endif()
endfunction()

function(expect_same_tree name actual reference)
  file(SHA256 "${actual}" actual_hash)
  file(SHA256 "${reference}" reference_hash)
  if(NOT actual_hash STREQUAL reference_hash)
    message(FATAL_ERROR "${name}: resumed tree differs from the uninterrupted one")
  endif()
endfunction()

foreach(mode default extended)
  set(mode_args)
  if(mode STREQUAL "extended")
    set(mode_args --extended)
  endif()
  set(ref "${output_dir}/${mode}-ref.tre")
  run_qraft("${mode} reference" TRUE "" ${mode_args} --threads 4 --out "${ref}")

  set(journal "${output_dir}/${mode}-full.journal")
  run_qraft("${mode} journaled" TRUE "" ${mode_args} --threads 4
            --out "${output_dir}/${mode}-full.tre"
            --checkpoint "${journal}" --checkpoint-interval 0)
  expect_same_tree("${mode} journaled" "${output_dir}/${mode}-full.tre" "${ref}")

  # A complete journal replays to the same tree without scoring.
  run_qraft("${mode} replay" TRUE "" ${mode_args} --threads 1
            --out "${output_dir}/${mode}-replay.tre"
            --checkpoint "${journal}" --resume)
  expect_same_tree("${mode} replay" "${output_dir}/${mode}-replay.tre" "${ref}")

  foreach(hook AFTER MID)
    foreach(n 1 2 4 8)
      set(case "${mode}-${hook}-${n}")
      set(journal "${output_dir}/${case}.journal")
      run_qraft("${case} cut" FALSE "QRAFT_CHECKPOINT_TEST_KILL_${hook}=${n}"
                ${mode_args} --threads 4 --out "${output_dir}/${case}-cut.tre"
                --checkpoint "${journal}" --checkpoint-interval 0)
      run_qraft("${case} resume" TRUE "" ${mode_args} --threads 2
                --out "${output_dir}/${case}.tre"
                --checkpoint "${journal}" --checkpoint-interval 0 --resume)
      expect_same_tree("${case}" "${output_dir}/${case}.tre" "${ref}")
    endforeach()
  endforeach()
endforeach()

# A long step does not hold back a due write. The run sleeps 60 s after its
# first decision (test hooks QRAFT_CHECKPOINT_TEST_STALL_AFTER and _MS), like
# a long step; with a 3 s interval the journal is written during the sleep,
# holding exactly that one decision, and the kill after that write ends the
# run there. The resumed run gives the uninterrupted tree.
set(journal "${output_dir}/stall.journal")
run_qraft("stall cut" FALSE
          "QRAFT_CHECKPOINT_TEST_STALL_AFTER=1;QRAFT_CHECKPOINT_TEST_STALL_MS=60000;QRAFT_CHECKPOINT_TEST_KILL_AFTER=1"
          --threads 4 --out "${output_dir}/stall-cut.tre"
          --checkpoint "${journal}" --checkpoint-interval 3)
file(STRINGS "${journal}" end_line REGEX "^END ")
if(NOT end_line MATCHES "^END 1 ")
  message(FATAL_ERROR
    "the journal write waited for the long step: first journal has '${end_line}'")
endif()
run_qraft("stall resume" TRUE "" --threads 2 --out "${output_dir}/stall.tre"
          --checkpoint "${journal}" --resume)
expect_same_tree("stall" "${output_dir}/stall.tre" "${output_dir}/default-ref.tre")

# A journal is not overwritten without --resume, and a journal of other
# options is refused.
run_qraft("existing journal" TRUE "" --threads 1
          --out "${output_dir}/guard-a.tre"
          --checkpoint "${output_dir}/guard.journal")
execute_process(
  COMMAND "${QRAFT_CLI}" --genes "${input}" --no-progress --quiet --threads 1
          --out "${output_dir}/guard-b.tre"
          --checkpoint "${output_dir}/guard.journal"
  RESULT_VARIABLE result ERROR_VARIABLE stderr)
if(result EQUAL 0 OR NOT stderr MATCHES "already exists")
  message(FATAL_ERROR "an existing journal was not protected:\n${stderr}")
endif()
execute_process(
  COMMAND "${QRAFT_CLI}" --genes "${input}" --no-progress --quiet --threads 1
          --no-nni --out "${output_dir}/guard-c.tre"
          --checkpoint "${output_dir}/guard.journal" --resume
  RESULT_VARIABLE result ERROR_VARIABLE stderr)
if(result EQUAL 0 OR NOT stderr MATCHES "options line differs")
  message(FATAL_ERROR "a journal of other options was accepted:\n${stderr}")
endif()

# A journal of other input is refused, even with the same taxa.
execute_process(
  COMMAND "${QRAFT_CLI}"
          --genes "${QRAFT_SOURCE_DIR}/tests/trees/astral2_200taxa_50genes_25missing.tre"
          --no-progress --quiet --threads 1 --out "${output_dir}/guard-d.tre"
          --checkpoint "${output_dir}/guard.journal" --resume
  RESULT_VARIABLE result ERROR_VARIABLE stderr)
if(result EQUAL 0 OR NOT stderr MATCHES "input line differs")
  message(FATAL_ERROR "a journal of other input was accepted:\n${stderr}")
endif()

# The status file may not be the journal.
execute_process(
  COMMAND "${QRAFT_CLI}" --genes "${input}" --no-progress --quiet --threads 1
          --out "${output_dir}/guard-e.tre"
          --checkpoint "${output_dir}/same.journal"
          --progress-file "${output_dir}/same.journal"
  RESULT_VARIABLE result ERROR_VARIABLE stderr)
if(result EQUAL 0 OR NOT stderr MATCHES "must differ")
  message(FATAL_ERROR "a status file on the journal was accepted:\n${stderr}")
endif()

# The status file is valid JSON and ends in state "done".
run_qraft("status file" TRUE "" --threads 2 --out "${output_dir}/status.tre"
          --progress-file "${output_dir}/status.json" --progress-interval 1)
file(READ "${output_dir}/status.json" status)
string(JSON state GET "${status}" state)
if(NOT state STREQUAL "done")
  message(FATAL_ERROR "status file ends in state ${state}")
endif()
