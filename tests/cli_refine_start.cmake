# --refine-start: one run gives the trees of the four runs that
# 'qraft help ablations' lists, for any worker count, with the run's options
# going to the run of the chain that uses them.

if(NOT DEFINED QRAFT_CLI OR NOT DEFINED QRAFT_SOURCE_DIR OR
   NOT DEFINED QRAFT_CLI_TEST_DIR)
  message(FATAL_ERROR "QRAFT_CLI, QRAFT_SOURCE_DIR, and QRAFT_CLI_TEST_DIR are required")
endif()

set(dir "${QRAFT_CLI_TEST_DIR}/refine-start")
file(REMOVE_RECURSE "${dir}")
file(MAKE_DIRECTORY "${dir}")

# Runs the CLI on ${genes} and writes ${dir}/${name}.tre.
function(run name)
  execute_process(
    COMMAND "${QRAFT_CLI}" --genes "${genes}" --out "${dir}/${name}.tre"
            --no-progress ${ARGN}
    RESULT_VARIABLE result
    OUTPUT_VARIABLE stdout
    ERROR_VARIABLE stderr)
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "${name} failed with ${result}\nstderr:\n${stderr}")
  endif()
endfunction()

function(same a b)
  file(SHA256 "${dir}/${a}.tre" hash_a)
  file(SHA256 "${dir}/${b}.tre" hash_b)
  if(NOT hash_a STREQUAL hash_b)
    message(FATAL_ERROR "${a}.tre and ${b}.tre differ")
  endif()
endfunction()

function(differ a b)
  file(SHA256 "${dir}/${a}.tre" hash_a)
  file(SHA256 "${dir}/${b}.tre" hash_b)
  if(hash_a STREQUAL hash_b)
    message(FATAL_ERROR "${a}.tre and ${b}.tre are the same")
  endif()
endfunction()

function(pinned name expected_hash)
  file(SHA256 "${dir}/${name}.tre" actual_hash)
  if(NOT actual_hash STREQUAL expected_hash)
    message(FATAL_ERROR
      "${name} hash changed: expected ${expected_hash}, got ${actual_hash}")
  endif()
endfunction()

# From the MRE consensus: the chain, then the flag. On these genes, which miss
# taxa, the refined start changes the output tree.
set(genes "${QRAFT_SOURCE_DIR}/tests/trees/astral2_200taxa_50genes_25missing.tre")
run(mre -t 4 --no-divide --no-nni)
run(refined -t 4 --starting-tree "${dir}/mre.tre" --no-divide --extended)
run(divided -t 4 --starting-tree "${dir}/refined.tre" --no-nni)
run(chain_default -t 4 --starting-tree "${dir}/divided.tre" --no-divide)
run(chain_extended -t 4 --starting-tree "${dir}/divided.tre" --no-divide
    --extended)
run(chain_support -t 4 --starting-tree "${dir}/divided.tre" --no-divide
    --support localpp)
run(divided_balanced -t 4 --starting-tree "${dir}/refined.tre" --no-nni
    --divide-score balanced)
run(default -t 4)
differ(refined mre)
differ(chain_default default)
pinned(chain_default
  4fa2eb77ff94f5966b17647d436ddea3eedf554ffc7075cef08b8a09ce0fc37e)
pinned(chain_extended
  163c54ebe7f47a7959d3075c0c1865deaa93ea3c693af219fcf9c4e07c79f001)
foreach(threads IN ITEMS 1 16)
  run(flag_default_${threads} -t ${threads} --refine-start
      --write-starting-tree "${dir}/flag_start_${threads}.tre")
  same(flag_default_${threads} chain_default)
  same(flag_start_${threads} refined)
  run(flag_extended_${threads} -t ${threads} --refine-start --extended)
  same(flag_extended_${threads} chain_extended)
  run(flag_no_nni_${threads} -t ${threads} --refine-start --no-nni)
  same(flag_no_nni_${threads} divided)
  run(flag_no_divide_${threads} -t ${threads} --refine-start --no-divide
      --no-nni)
  same(flag_no_divide_${threads} refined)
endforeach()
run(flag_support -t 4 --refine-start --support localpp)
same(flag_support chain_support)
run(flag_balanced -t 4 --refine-start --no-nni --divide-score balanced)
same(flag_balanced divided_balanced)

# From a given starting tree with planted errors, and with it as the guide.
set(genes "${QRAFT_SOURCE_DIR}/tests/trees/astral2_200taxa_1000genes.tre")
set(planted "${QRAFT_SOURCE_DIR}/tests/trees/astral2_200taxa_1000genes_planted.tre")
run(planted_start -t 4 --starting-tree "${planted}" --no-divide --no-nni)
run(planted_refined -t 4 --starting-tree "${dir}/planted_start.tre"
    --no-divide --extended)
run(planted_divided -t 4 --starting-tree "${dir}/planted_refined.tre" --no-nni)
run(planted_chain -t 4 --starting-tree "${dir}/planted_divided.tre" --no-divide)
differ(planted_refined planted_start)
run(planted_flag -t 4 --refine-start --starting-tree "${planted}")
same(planted_flag planted_chain)
run(guided_start -t 4 --guide-tree "${planted}" --no-divide --no-nni)
run(guided_refined -t 4 --starting-tree "${dir}/guided_start.tre" --no-divide
    --extended)
run(guided_flag -t 4 --refine-start --guide-tree "${planted}" --no-divide
    --no-nni)
same(guided_flag guided_refined)

# The statistics and the run summary name the setting.
set(genes "${QRAFT_SOURCE_DIR}/tests/trees/astral2_200taxa_50genes.tre")
run(stats -t 2 --refine-start --no-nni -v --stats-json "${dir}/stats.json")
file(READ "${dir}/stats.json" stats)
string(JSON refine_start GET "${stats}" refine_start)
if(NOT refine_start STREQUAL "ON" AND NOT refine_start STREQUAL "true")
  message(FATAL_ERROR "refine_start missing from the statistics:\n${stats}")
endif()
execute_process(
  COMMAND "${QRAFT_CLI}" --genes "${genes}" --out "${dir}/summary.tre"
          --no-progress -v -t 2 --refine-start --no-nni
  RESULT_VARIABLE result
  OUTPUT_VARIABLE stdout
  ERROR_VARIABLE stderr)
if(NOT result EQUAL 0 OR NOT stderr MATCHES "custom settings, refined starting tree")
  message(FATAL_ERROR "the run summary does not name --refine-start:\n${stderr}")
endif()
