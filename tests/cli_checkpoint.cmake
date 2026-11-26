# Checks the decision journal (--checkpoint, --resume): a journaled run gives
# the same tree as a plain run, and a run killed right after, or in the middle
# of, a journal write and then resumed gives that tree too, at another thread
# count. The kills use the test hooks QRAFT_CHECKPOINT_TEST_KILL_AFTER and
# QRAFT_CHECKPOINT_TEST_KILL_MID, which end the process like kill -9. A cut
# NNI invocation runs again from its start, also inside the SPR refinement,
# and a due journal write is not held back by a long step.

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

# Cuts inside NNI and SPR. With a kill hook and --checkpoint-interval 0 each
# decision waits for the journal write that holds it, so
# QRAFT_CHECKPOINT_TEST_KILL_AFTER=n cuts right after the n-th record. From a
# starting tree with --no-divide the journal holds only NNI and SPR records,
# in a fixed order. A cut NNI invocation runs again from its start: the
# journal written with the resumed run's first record holds the records
# before that invocation and then its first record again. An invocation cut
# after its end is replayed without scoring. Every resumed run gives the
# uninterrupted tree and journal.

# The journal's records, without its header and END line.
function(journal_records out path)
  file(STRINGS "${path}" records REGEX "^[DGMNPQS]( |$)")
  set(${out} "${records}" PARENT_SCOPE)
endfunction()

function(expect_records name path expected)
  journal_records(actual "${path}")
  if(NOT actual STREQUAL expected)
    list(LENGTH actual actual_count)
    list(LENGTH expected expected_count)
    message(FATAL_ERROR
      "${name}: the journal holds ${actual_count} records, not the "
      "${expected_count} expected, or other ones")
  endif()
endfunction()

function(run_start name expect_success env genes start)
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E env ${env}
            "${QRAFT_CLI}" --genes "${genes}" --starting-tree "${start}"
            --no-divide --no-progress --quiet ${ARGN}
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

# run_cuts(TAG GENES START ARGS <options> CUTS <n...> AT <kind inv>...):
# AT names the kind and invocation of each cut's record ("N 1", "M 0").
function(run_cuts tag genes start)
  cmake_parse_arguments(arg "" "" "ARGS;CUTS;AT" ${ARGN})
  set(ref "${output_dir}/${tag}-ref.tre")
  run_start("${tag} reference" TRUE "" "${genes}" "${start}" ${arg_ARGS}
            --threads 4 --out "${ref}")
  set(full "${output_dir}/${tag}-full.journal")
  run_start("${tag} journaled" TRUE "" "${genes}" "${start}" ${arg_ARGS}
            --threads 4 --out "${output_dir}/${tag}-full.tre"
            --checkpoint "${full}" --checkpoint-interval 0)
  expect_same_tree("${tag} journaled" "${output_dir}/${tag}-full.tre" "${ref}")
  journal_records(records "${full}")
  list(LENGTH arg_CUTS cut_count)
  math(EXPR last "${cut_count} - 1")
  foreach(i RANGE ${last})
    list(GET arg_CUTS ${i} n)
    list(GET arg_AT ${i} at)
    math(EXPR index "${n} - 1")
    list(GET records ${index} record)
    if(NOT record MATCHES "^${at}( |$)")
      message(FATAL_ERROR "${tag}: record ${n} is '${record}', not ${at}")
    endif()
    list(SUBLIST records 0 ${n} before)
    # The records the resumed run keeps: all before the cut NNI invocation.
    set(kept ${n})
    if(record MATCHES "^N ([0-9]+) ")
      set(invocation "${CMAKE_MATCH_1}")
      set(kept 0)
      foreach(other IN LISTS records)
        if(other MATCHES "^N ${invocation} ")
          break()
        endif()
        math(EXPR kept "${kept} + 1")
      endforeach()
    endif()
    math(EXPR kept_and_next "${kept} + 1")
    list(SUBLIST records 0 ${kept_and_next} after_next)

    set(case "${tag}-${n}")
    set(journal "${output_dir}/${case}.journal")
    run_start("${case} cut" FALSE "QRAFT_CHECKPOINT_TEST_KILL_AFTER=${n}"
              "${genes}" "${start}" ${arg_ARGS} --threads 4
              --out "${output_dir}/${case}-cut.tre"
              --checkpoint "${journal}" --checkpoint-interval 0)
    expect_records("${case} cut" "${journal}" "${before}")
    run_start("${case} resume" TRUE "" "${genes}" "${start}" ${arg_ARGS}
              --threads 2 --out "${output_dir}/${case}.tre"
              --checkpoint "${journal}" --checkpoint-interval 0 --resume)
    expect_same_tree("${case}" "${output_dir}/${case}.tre" "${ref}")
    expect_records("${case} resume" "${journal}" "${records}")

    set(journal "${output_dir}/${case}-twice.journal")
    run_start("${case} cut once" FALSE "QRAFT_CHECKPOINT_TEST_KILL_AFTER=${n}"
              "${genes}" "${start}" ${arg_ARGS} --threads 4
              --out "${output_dir}/${case}-cut1.tre"
              --checkpoint "${journal}" --checkpoint-interval 0)
    run_start("${case} cut twice" FALSE "QRAFT_CHECKPOINT_TEST_KILL_AFTER=1"
              "${genes}" "${start}" ${arg_ARGS} --threads 2
              --out "${output_dir}/${case}-cut2.tre"
              --checkpoint "${journal}" --checkpoint-interval 0 --resume)
    expect_records("${case} cut twice" "${journal}" "${after_next}")
    run_start("${case} resume twice" TRUE "" "${genes}" "${start}" ${arg_ARGS}
              --threads 2 --out "${output_dir}/${case}-twice.tre"
              --checkpoint "${journal}" --checkpoint-interval 0 --resume)
    expect_same_tree("${case} twice" "${output_dir}/${case}-twice.tre" "${ref}")
    expect_records("${case} resume twice" "${journal}" "${records}")
  endforeach()
endfunction()

# The planted 200-taxon tree, 1,000 genes: default NNI makes 5 single moves
# (records 1-5) and --nni-batch a few batches.
set(planted "${QRAFT_SOURCE_DIR}/tests/trees/astral2_200taxa_1000genes_planted.tre")
set(planted_genes "${QRAFT_SOURCE_DIR}/tests/trees/astral2_200taxa_1000genes.tre")
run_cuts(nni-default "${planted_genes}" "${planted}"
         CUTS 1 2 4 AT "N 0" "N 0" "N 0")
run_cuts(nni-batch "${planted_genes}" "${planted}" ARGS --nni-batch
         CUTS 1 2 AT "N 0" "N 0")

# A random starting tree, 50 genes. With --spr-refine, NNI 0 makes 20 moves
# (records 1-20, its end 21), SPR round 0 (22-32), nested NNI 1 (33-52, end
# 53), round 1 (54-70), nested NNI 2 (71-90, end 91). With --nni-batch too,
# NNI 0 and its end (1-2), round 0 (3-15), nested NNI 1 one batch (16, end
# 17). Extended mode: NNI 0 makes 55 batches (1-55, end 56), and SPR finds
# nothing to move.
set(random_start
    "${QRAFT_SOURCE_DIR}/tests/trees/astral2_200taxa_50genes_random_start.tre")
run_cuts(spr "${input}" "${random_start}" ARGS --spr-refine
         CUTS 25 40 53 60 80 AT "P 0" "N 1" "M 1" "P 0" "N 2")
run_cuts(spr-batch "${input}" "${random_start}" ARGS --nni-batch --spr-refine
         CUTS 16 AT "N 1")
run_cuts(random-extended "${input}" "${random_start}" ARGS --extended
         CUTS 30 56 AT "N 0" "M 0")

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
