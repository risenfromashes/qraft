# The SPR refinement through the CLI: a planted misplaced clade is moved home
# on a starting tree, with the same output at 1 and 4 workers and under either
# regraft score; invalid parameters are refused; --extended sets the four
# options; a tree without eligible clades reads nothing; the divide path runs
# the stage too.

if(NOT DEFINED QRAFT_CLI OR NOT DEFINED QRAFT_CLI_TEST_DIR)
  message(FATAL_ERROR "QRAFT_CLI and QRAFT_CLI_TEST_DIR are required")
endif()

set(dir "${QRAFT_CLI_TEST_DIR}/spr-cli")
file(MAKE_DIRECTORY "${dir}")
file(WRITE "${dir}/genes.tre"
  "(((A,B),(C,D)),((E,F),(G,H)));\n(((A,B),(C,D)),((E,F),(G,H)));\n")
file(WRITE "${dir}/start.tre" "((C,D),(((A,B),(E,F)),(G,H)));\n")
foreach(workers 1 4)
  foreach(mode raw normalized)
    set(stem "${dir}/${workers}-${mode}")
    execute_process(COMMAND "${QRAFT_CLI}" -i "${dir}/genes.tre"
      --starting-tree "${dir}/start.tre" --no-divide --no-nni
      --spr-refine --spr-min-clade 2
      --regraft-score ${mode} -t ${workers} --no-progress
      --stats-json "${stem}.json" -o "${stem}.tre" RESULT_VARIABLE result)
    if(NOT result EQUAL 0)
      message(FATAL_ERROR "SPR CLI failed: ${result}")
    endif()
    file(READ "${stem}.json" stats)
    string(JSON moves GET "${stats}" spr_moves)
    string(JSON gain GET "${stats}" spr_gain)
    string(JSON reads GET "${stats}" spr_reads)
    if(moves LESS 1 OR gain LESS_EQUAL 0 OR reads LESS 1)
      message(FATAL_ERROR "SPR refinement applied no improvement:\n${stats}")
    endif()
    file(READ "${stem}.tre" tree)
    if(DEFINED reference AND NOT tree STREQUAL reference)
      message(FATAL_ERROR "SPR depends on the worker count or the regraft score")
    endif()
    set(reference "${tree}")
  endforeach()
endforeach()
foreach(flag --spr-reads --spr-min-clade --spr-candidates --spr-rounds)
  foreach(value 0 -1 nope)
    execute_process(COMMAND "${QRAFT_CLI}" -i "${dir}/genes.tre" ${flag} ${value}
      RESULT_VARIABLE result OUTPUT_QUIET ERROR_QUIET)
    if(result EQUAL 0)
      message(FATAL_ERROR "SPR accepted invalid ${flag} ${value}")
    endif()
  endforeach()
endforeach()
execute_process(COMMAND "${QRAFT_CLI}" -i "${dir}/genes.tre" --spr-min-clade 1
  RESULT_VARIABLE result OUTPUT_QUIET ERROR_QUIET)
if(result EQUAL 0)
  message(FATAL_ERROR "SPR accepted a minimum clade of 1")
endif()
# --extended is unbounded batched NNI and the SPR refinement.
execute_process(COMMAND "${QRAFT_CLI}" -i "${dir}/genes.tre" --extended
  --spr-min-clade 2 -t 1 --no-progress --stats-json "${dir}/extended.json"
  -o "${dir}/extended.tre" RESULT_VARIABLE result)
if(NOT result EQUAL 0)
  message(FATAL_ERROR "extended mode failed: ${result}")
endif()
file(READ "${dir}/extended.json" stats)
foreach(key nni_batch nni_move_cap nni_support_gate spr_refine spr_rounds)
  string(JSON ${key} GET "${stats}" ${key})
endforeach()
if(NOT nni_batch OR nni_move_cap OR nni_support_gate OR NOT spr_refine OR
   spr_rounds LESS 1)
  message(FATAL_ERROR "--extended did not set its options:\n${stats}")
endif()
# At the default minimum clade no clade of eight taxa is eligible: the round
# reads nothing, moves nothing and ends the refinement.
execute_process(COMMAND "${QRAFT_CLI}" -i "${dir}/genes.tre"
  --starting-tree "${dir}/start.tre" --no-divide --no-nni --spr-refine
  -t 1 --no-progress --stats-json "${dir}/ineligible.json"
  -o "${dir}/ineligible.tre" RESULT_VARIABLE result)
if(NOT result EQUAL 0)
  message(FATAL_ERROR "SPR without eligible clades failed: ${result}")
endif()
file(READ "${dir}/ineligible.json" stats)
foreach(key spr_rounds spr_reads spr_candidates_scored spr_moves)
  string(JSON ${key} GET "${stats}" ${key})
endforeach()
if(NOT spr_rounds EQUAL 1 OR NOT spr_reads EQUAL 0 OR
   NOT spr_candidates_scored EQUAL 0 OR NOT spr_moves EQUAL 0)
  message(FATAL_ERROR "SPR read or moved without an eligible clade:\n${stats}")
endif()
# The divide path runs the stage too, even without NNI.
execute_process(COMMAND "${QRAFT_CLI}" -i "${dir}/genes.tre" --no-nni
  --spr-refine --spr-min-clade 2 -t 1 --no-progress
  --stats-json "${dir}/divide.json" -o "${dir}/divide.tre"
  RESULT_VARIABLE result)
if(NOT result EQUAL 0)
  message(FATAL_ERROR "SPR divide path failed: ${result}")
endif()
file(READ "${dir}/divide.json" stats)
string(JSON rounds GET "${stats}" spr_rounds)
if(rounds LESS 1)
  message(FATAL_ERROR "SPR divide path did not run the stage")
endif()
