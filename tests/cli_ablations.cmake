# Pins the output hash of the default pipeline and of each ablation, and
# checks that the worker count does not change the output.

if(NOT DEFINED QRAFT_CLI OR NOT DEFINED QRAFT_SOURCE_DIR OR
   NOT DEFINED QRAFT_CLI_TEST_DIR)
  message(FATAL_ERROR "QRAFT_CLI, QRAFT_SOURCE_DIR, and QRAFT_CLI_TEST_DIR are required")
endif()

set(output_dir "${QRAFT_CLI_TEST_DIR}/supported-mode-hashes")
file(MAKE_DIRECTORY "${output_dir}")
set(threads 4)

# Runs the CLI on ${input} with ${threads} workers and checks the output hash.
function(assert_mode_hash name expected_hash)
  set(output "${output_dir}/${name}.tre")
  execute_process(
    COMMAND "${QRAFT_CLI}"
            --genes "${input}"
            --out "${output}"
            --threads ${threads}
            --no-progress
            ${ARGN}
    RESULT_VARIABLE result
    OUTPUT_VARIABLE stdout
    ERROR_VARIABLE stderr)
  if(NOT result EQUAL 0)
    message(FATAL_ERROR
      "${name} failed with ${result}\nstdout:\n${stdout}\nstderr:\n${stderr}")
  endif()
  file(SHA256 "${output}" actual_hash)
  if(NOT actual_hash STREQUAL expected_hash)
    message(FATAL_ERROR
      "${name} hash changed: expected ${expected_hash}, got ${actual_hash}")
  endif()
endfunction()

# 201 taxa, 1000 genes.
set(input "${QRAFT_SOURCE_DIR}/tests/trees/astral2_200taxa_1000genes.tre")
assert_mode_hash(qraft
  b6c6d8b0d05c75d3d0d15f973c3dae01535f80b97a5c55c170c55e9db427b30e)
assert_mode_hash(no_bipartition_test
  b6c6d8b0d05c75d3d0d15f973c3dae01535f80b97a5c55c170c55e9db427b30e
  --no-bipartition-test)
# A crossover below the gene count makes the threshold positive, so the test
# rejects bipartitions; more proposals and anchor pairs change the outcome.
assert_mode_hash(bipartition_test_rejecting
  db32c17b3ac0a2d4fda143f39a13cc7477b374c7ad950c3f3da84b7863bff494
  --test-crossover 100)
assert_mode_hash(bipartition_test_retries
  47a15d345c03d4a50948035b7c0493825ac2651d6c8ab95ae6271e59c14b9fcd
  --test-crossover 100 --test-proposals 4 --test-anchor-pairs 4)
assert_mode_hash(balanced_divide
  728a17b9fb5a78cb4eaabc585e9519b92e02fd0ec5a02c6e40e880335f1d79d9
  --divide-score balanced)
assert_mode_hash(starting_tree_nni
  9c4e7f86852aedac7b29152c739d4688865ddc5a42a07b37498d094632d0672c
  --no-divide)
assert_mode_hash(starting_tree
  5add8f4e027f9b56768a54c625d6937cbf833439588276a7cf8e849b889cc2ba
  --no-divide --no-nni)
# Disabling the bipartition test must not build its distance index.
assert_mode_hash(balanced_divide_without_test
  728a17b9fb5a78cb4eaabc585e9519b92e02fd0ec5a02c6e40e880335f1d79d9
  --divide-score balanced --no-bipartition-test
  --stats-json "${output_dir}/without-test.json")
file(READ "${output_dir}/without-test.json" without_test)
string(JSON test_available GET "${without_test}" bipartition_test_available)
string(JSON test_count GET "${without_test}" bipartition_tests)
if(test_available OR NOT test_count EQUAL 0)
  message(FATAL_ERROR "the bipartition-test ablation still built its index")
endif()

# 200 taxa, 50 genes: NNI moves, and each ablation changes the tree.
set(input "${QRAFT_SOURCE_DIR}/tests/trees/astral2_200taxa_50genes.tre")
assert_mode_hash(qraft_50_genes
  f0b8a1d64ad08890c060eb19026f107942de6532adde6e02e03e6a4c2c5475c7)
assert_mode_hash(normalized_regraft_50_genes
  9d4d4e5c82609f8ab376a43b626bfb89d5d4aad824a9ceb99bf13fdc2a49293c
  --regraft-score normalized)
assert_mode_hash(no_nni_50_genes
  0afc0b574b5366d41d9260564fdbf9c4f53e354758ad9373c00e5207ad1580b3
  --no-nni)
assert_mode_hash(balanced_divide_50_genes
  5792f1329454ca23f318c5d6e05fa6c91abd4b41ce8ec156830f19910a087fe8
  --divide-score balanced)
assert_mode_hash(starting_tree_50_genes
  eb5f099600dffcfa6b74228f6e2ba68965dbfaecc4ae64684675eeb3321b8561
  --no-divide)
assert_mode_hash(no_nni_support_gate_50_genes
  2544a2b4729a0e80c34d017eae4573014569cc79a078e6da0599e0bb44ce5378
  --no-nni-support-gate)
# The batch skips the moves the support gate rejects instead of stopping.
assert_mode_hash(nni_batch_50_genes
  48f250be9de274dd7ca5140193a52f7589a80c7b117489f8d08cacfff3614d66
  --nni-batch)
assert_mode_hash(nni_batch_ungated_uncapped_50_genes
  2544a2b4729a0e80c34d017eae4573014569cc79a078e6da0599e0bb44ce5378
  --nni-batch --no-nni-cap --no-nni-support-gate)

# 200 taxa, 1000 genes, NNI on the starting tree: without the support gate
# NNI reaches the cap of 20 moves, and makes 27 without it; one batch pass
# offers more moves than the cap leaves, so the batch keeps its largest gains.
set(input "${QRAFT_SOURCE_DIR}/tests/trees/astral2_200taxa_1000genes.tre")
assert_mode_hash(starting_tree_ungated
  fc00dd4cba976875c363a7c5b1ac8fc0292dd262fcbe9ad3f1b631b2d9d893af
  --no-divide --no-nni-support-gate)
assert_mode_hash(starting_tree_ungated_uncapped
  d63ed76a37589637fb56ada52ed0e6d222bcc0667ab7257dfaab20f7f426172a
  --no-divide --no-nni-support-gate --no-nni-cap)
assert_mode_hash(starting_tree_batch
  d94db46663d4217dfa56679a06bb327c7b3575d41add6d9e891fe2f427a5b1b6
  --no-divide --nni-batch)
assert_mode_hash(starting_tree_batch_ungated
  5e61174d1665546c298bf6111595d2ae7e116c56c9675ae73e88320469cb81aa
  --no-divide --nni-batch --no-nni-support-gate)
assert_mode_hash(starting_tree_batch_ungated_uncapped
  d63ed76a37589637fb56ada52ed0e6d222bcc0667ab7257dfaab20f7f426172a
  --no-divide --nni-batch --no-nni-support-gate --no-nni-cap
  --stats-json "${output_dir}/batch.json")
# The batch statistics, and no support cache without the gate.
file(READ "${output_dir}/batch.json" batch_stats)
foreach(key nni_batch nni_move_cap nni_support_gate nni_moves_accepted
        nni_batch_passes nni_gain nni_support_stops nni_split_support_ns)
  string(JSON ${key} GET "${batch_stats}" ${key})
endforeach()
if(NOT nni_batch OR nni_move_cap OR nni_support_gate OR
   NOT nni_moves_accepted EQUAL 27 OR NOT nni_batch_passes GREATER 0 OR
   NOT nni_gain EQUAL 31797877 OR NOT nni_support_stops EQUAL 0 OR
   NOT nni_split_support_ns EQUAL 0)
  message(FATAL_ERROR "unexpected batch NNI statistics:\n${batch_stats}")
endif()
assert_mode_hash(starting_tree_nni_stats
  9c4e7f86852aedac7b29152c739d4688865ddc5a42a07b37498d094632d0672c
  --no-divide --stats-json "${output_dir}/sequential.json")
file(READ "${output_dir}/sequential.json" sequential_stats)
string(JSON sequential_passes GET "${sequential_stats}" nni_batch_passes)
string(JSON sequential_gain GET "${sequential_stats}" nni_gain)
if(NOT sequential_passes EQUAL 0 OR NOT sequential_gain EQUAL 23839542)
  message(FATAL_ERROR "unexpected sequential NNI statistics:\n${sequential_stats}")
endif()

# Extended mode. On the simulated inputs the SPR pass finds no clade with a
# positive skeleton gain, so --extended equals unbounded batched NNI.
set(input "${QRAFT_SOURCE_DIR}/tests/trees/astral2_200taxa_50genes.tre")
assert_mode_hash(extended_50_genes
  2544a2b4729a0e80c34d017eae4573014569cc79a078e6da0599e0bb44ce5378
  --extended)
# A starting tree with a misplaced 21-taxon clade: the SPR refinement alone
# (at 1 and 16 workers, below) applies three exact re-joins in its first
# round and none in its second; under --extended, unbounded NNI repairs the
# tree first and the SPR pass finds no candidate.
set(input "${QRAFT_SOURCE_DIR}/tests/trees/astral2_200taxa_1000genes.tre")
set(planted "${QRAFT_SOURCE_DIR}/tests/trees/astral2_200taxa_1000genes_planted.tre")
assert_mode_hash(extended_planted
  fc1508423097ef2c9645848636e54a1391621bacf34ef4e45b44f23da7ecbe3f
  --starting-tree "${planted}" --no-divide --extended --spr-min-clade 10)

# 200 taxa, 50 genes that each miss a quarter of their taxa: the bipartition
# test runs on the incomplete genes with J on shared genes. With 50 genes its
# threshold is lenient and it accepts every bipartition; with a crossover of 20
# it rejects most of them, and with a crossover of 50 it both accepts and
# rejects.
set(input "${QRAFT_SOURCE_DIR}/tests/trees/astral2_200taxa_50genes_25missing.tre")
assert_mode_hash(incomplete_genes
  a79188c47d889c560fe5535cebe519271605c02be096247a54e49105ff1866eb)
assert_mode_hash(incomplete_genes_rejecting
  4eb74beafe3668a44b3e971a0fda801ca8df53cec05f693e418398470ec95ce0
  --test-crossover 20 --stats-json "${output_dir}/incomplete.json")
file(READ "${output_dir}/incomplete.json" incomplete_stats)
foreach(key bipartition_test_available bipartition_test_complete_genes
        bipartition_test_incomplete bipartition_test_rejects)
  string(JSON ${key} GET "${incomplete_stats}" ${key})
endforeach()
if(NOT bipartition_test_available OR bipartition_test_complete_genes OR
   NOT bipartition_test_incomplete GREATER 0 OR
   NOT bipartition_test_rejects GREATER 0)
  message(FATAL_ERROR "the bipartition test did not run on incomplete genes:\n${incomplete_stats}")
endif()
assert_mode_hash(incomplete_genes_mixed
  fff47bcc7eabb7e4cb96fedab23305928fa0deef6a77ac5203ee10960fcbc1ce
  --test-crossover 50 --stats-json "${output_dir}/incomplete-mixed.json")
file(READ "${output_dir}/incomplete-mixed.json" mixed_stats)
string(JSON mixed_accepts GET "${mixed_stats}" bipartition_test_accepts)
string(JSON mixed_rejects GET "${mixed_stats}" bipartition_test_rejects)
if(NOT mixed_accepts GREATER 0 OR NOT mixed_rejects GREATER 0)
  message(FATAL_ERROR "the incomplete-gene test did not both accept and reject:\n${mixed_stats}")
endif()
foreach(threads IN ITEMS 1 16)
  assert_mode_hash(incomplete_genes_rejecting_threads_${threads}
    4eb74beafe3668a44b3e971a0fda801ca8df53cec05f693e418398470ec95ce0
    --test-crossover 20)
  assert_mode_hash(incomplete_genes_mixed_threads_${threads}
    fff47bcc7eabb7e4cb96fedab23305928fa0deef6a77ac5203ee10960fcbc1ce
    --test-crossover 50)
endforeach()
set(threads 4)

# The worker count does not change the output.
set(input "${QRAFT_SOURCE_DIR}/tests/trees/astral2_200taxa_50genes.tre")
foreach(threads IN ITEMS 1 16)
  assert_mode_hash(qraft_50_genes_threads_${threads}
    f0b8a1d64ad08890c060eb19026f107942de6532adde6e02e03e6a4c2c5475c7)
  assert_mode_hash(nni_batch_50_genes_threads_${threads}
    48f250be9de274dd7ca5140193a52f7589a80c7b117489f8d08cacfff3614d66
    --nni-batch)
endforeach()
set(threads 4)
set(input "${QRAFT_SOURCE_DIR}/tests/trees/astral2_200taxa_1000genes.tre")
foreach(threads IN ITEMS 1 16)
  assert_mode_hash(starting_tree_nni_threads_${threads}
    9c4e7f86852aedac7b29152c739d4688865ddc5a42a07b37498d094632d0672c
    --no-divide)
  assert_mode_hash(starting_tree_batch_ungated_uncapped_threads_${threads}
    d63ed76a37589637fb56ada52ed0e6d222bcc0667ab7257dfaab20f7f426172a
    --no-divide --nni-batch --no-nni-support-gate --no-nni-cap)
  assert_mode_hash(spr_refine_planted_threads_${threads}
    7df6f08820c547392c6ba583e2d7b29ad3197a390ded66f25241ed7379ab8da2
    --starting-tree "${planted}" --no-divide --no-nni --spr-refine
    --spr-min-clade 10 --stats-json "${output_dir}/spr_${threads}.json")
  file(READ "${output_dir}/spr_${threads}.json" spr_stats)
  foreach(key spr_rounds spr_candidates_scored spr_moves spr_gain)
    string(JSON ${key} GET "${spr_stats}" ${key})
  endforeach()
  if(NOT spr_rounds EQUAL 2 OR NOT spr_candidates_scored EQUAL 9 OR
     NOT spr_moves EQUAL 3 OR NOT spr_gain EQUAL 3778327051)
    message(FATAL_ERROR "unexpected SPR refinement statistics:\n${spr_stats}")
  endif()
endforeach()
set(threads 4)

execute_process(
  COMMAND "${QRAFT_CLI}"
          --genes "${input}"
          --out "${output_dir}/unknown-option.tre"
          --no-progress
          --no-such-option
  RESULT_VARIABLE result
  OUTPUT_VARIABLE stdout
  ERROR_VARIABLE stderr)
if(NOT result EQUAL 2 OR NOT stderr MATCHES "unknown argument: --no-such-option")
  message(FATAL_ERROR "an unknown option was not a usage error\n"
    "result=${result}\nstdout:\n${stdout}\nstderr:\n${stderr}")
endif()
