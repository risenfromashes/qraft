# Runs the CLI with a 256 KiB stack, at one and at four threads, on inputs
# that nest deeply. Caterpillar gene trees nest as deep as they have taxa in
# the Newick parser, and their MRE starting tree is one polytomy that the
# divide peels one arm per level; so is the MRE of stars of cherries whose
# pairings disagree. The expected hashes are those of the trees the earlier
# recursive implementation wrote with a large stack; -DQRAFT_STACK_KIB=N runs
# with another stack size, for regenerating them.

if(NOT DEFINED QRAFT_CLI OR NOT DEFINED QRAFT_CLI_TEST_DIR)
  message(FATAL_ERROR "QRAFT_CLI and QRAFT_CLI_TEST_DIR are required")
endif()

set(work_dir "${QRAFT_CLI_TEST_DIR}/deep-inputs")
file(MAKE_DIRECTORY "${work_dir}")
if(DEFINED QRAFT_STACK_KIB)
  set(stack_kib ${QRAFT_STACK_KIB})
else()
  set(stack_kib 256)
endif()

# Writes `genes` stars of cherries on `taxa` taxa (an even number): even genes
# pair t(2j) with t(2j+1), odd genes t(2j+1) with t(2j+2), wrapping around,
# so every cherry is in half of the genes.
function(write_cherry_stars path taxa genes)
  math(EXPR last_pair "${taxa} / 2 - 1")
  math(EXPR last_gene "${genes} - 1")
  set(text "")
  foreach(gene RANGE ${last_gene})
    math(EXPR shift "${gene} % 2")
    set(cherries "")
    foreach(pair RANGE ${last_pair})
      math(EXPR first "(2 * ${pair} + ${shift}) % ${taxa}")
      math(EXPR second "(2 * ${pair} + ${shift} + 1) % ${taxa}")
      list(APPEND cherries "(t${first},t${second})")
    endforeach()
    list(JOIN cherries "," joined)
    string(APPEND text "(${joined});\n")
  endforeach()
  file(WRITE "${path}" "${text}")
endfunction()

# Writes two caterpillars on `taxa` taxa: one in taxon order, one in the
# order i * 7919 mod taxa, which is a permutation when taxa is coprime to
# 7919.
function(write_caterpillars path taxa)
  math(EXPR last "${taxa} - 1")
  math(EXPR opening "${taxa} - 1")
  string(REPEAT "(" ${opening} prefix)
  set(text "")
  foreach(stride 1 7919)
    set(body "")
    foreach(i RANGE 1 ${last})
      math(EXPR taxon "(${i} * ${stride}) % ${taxa}")
      string(APPEND body ",t${taxon})")
    endforeach()
    string(APPEND text "${prefix}t0${body};\n")
  endforeach()
  file(WRITE "${path}" "${text}")
endfunction()

# Runs the CLI on ${input} with a small stack at 1 and 4 threads and checks
# the output hash.
function(assert_small_stack_hash name input expected_hash)
  foreach(threads 1 4)
    set(output "${work_dir}/${name}-t${threads}.tre")
    execute_process(
      COMMAND sh -c "ulimit -s ${stack_kib} && exec \"$0\" \"$@\""
              "${QRAFT_CLI}"
              --genes "${input}"
              --out "${output}"
              --threads ${threads}
              --no-progress
      RESULT_VARIABLE result
      OUTPUT_VARIABLE stdout
      ERROR_VARIABLE stderr)
    if(NOT result EQUAL 0)
      message(FATAL_ERROR
        "${name} at ${threads} threads with a ${stack_kib} KiB stack failed "
        "with ${result}\nstdout:\n${stdout}\nstderr:\n${stderr}")
    endif()
    file(SHA256 "${output}" actual_hash)
    if(NOT actual_hash STREQUAL expected_hash)
      message(FATAL_ERROR
        "${name} at ${threads} threads: expected ${expected_hash}, "
        "got ${actual_hash}")
    endif()
  endforeach()
endfunction()

set(stars "${work_dir}/cherry-stars-600x8.tre")
set(caterpillars "${work_dir}/caterpillars-10000x2.tre")
write_cherry_stars("${stars}" 600 8)
write_caterpillars("${caterpillars}" 10000)
assert_small_stack_hash(cherry_stars "${stars}"
  c7233c90f38779161b9af56d9125754f1361bf311d10e632be2a7e48ae91335d)
assert_small_stack_hash(caterpillars "${caterpillars}"
  0beb649b42ea9044501fb18bb5b230ea9fb29829f1b0ac005d08f763de5bcf3f)
