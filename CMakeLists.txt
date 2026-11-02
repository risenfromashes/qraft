cmake_minimum_required(VERSION 3.20)

project(qraft VERSION 1.0.0 LANGUAGES CXX)

option(QRAFT_BUILD_TESTS "Build QRAFT tests" ON)
option(QRAFT_REGISTER_LONG_TESTS
  "Register exhaustive long-running proof shards" OFF)
option(QRAFT_NATIVE_OPTIMIZATIONS
  "Build release binaries with machine-local native CPU tuning" OFF)
option(QRAFT_ENABLE_IPO "Enable interprocedural optimization for release builds" ON)
option(QRAFT_USE_JEMALLOC
  "Use jemalloc; required by the default static executable profile" ON)
option(QRAFT_STATIC_EXECUTABLE
  "Build qraft as a fully static executable (Linux, or Windows with MinGW-w64)" ON)

set(QRAFT_NATIVE_CPU_FLAG_PATTERN
  "(^|[^A-Za-z0-9_])-(march|mcpu|mtune)=native($|[^A-Za-z0-9_])")
string(CONCAT QRAFT_CONFIGURED_CXX_AND_LINK_FLAGS
  "${CMAKE_CXX_COMPILER_ARG1} "
  "${CMAKE_CXX_FLAGS} "
  "${CMAKE_CXX_FLAGS_DEBUG} "
  "${CMAKE_CXX_FLAGS_RELEASE} "
  "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} "
  "${CMAKE_CXX_FLAGS_MINSIZEREL} "
  "${CMAKE_EXE_LINKER_FLAGS} "
  "${CMAKE_EXE_LINKER_FLAGS_DEBUG} "
  "${CMAKE_EXE_LINKER_FLAGS_RELEASE} "
  "${CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO} "
  "${CMAKE_EXE_LINKER_FLAGS_MINSIZEREL}")

if(QRAFT_STATIC_EXECUTABLE)
  if(QRAFT_NATIVE_OPTIMIZATIONS)
    message(FATAL_ERROR
      "QRAFT_STATIC_EXECUTABLE requires QRAFT_NATIVE_OPTIMIZATIONS=OFF; "
      "a machine-local instruction set is not a portable artifact")
  endif()
  if(BUILD_SHARED_LIBS)
    message(FATAL_ERROR
      "QRAFT_STATIC_EXECUTABLE requires BUILD_SHARED_LIBS=OFF")
  endif()
  if(NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$" OR
     NOT (CMAKE_SYSTEM_NAME STREQUAL "Linux" OR MINGW))
    message(FATAL_ERROR
      "QRAFT_STATIC_EXECUTABLE requires Linux, or Windows with MinGW-w64, and "
      "GNU or Clang")
  endif()
  if(QRAFT_CONFIGURED_CXX_AND_LINK_FLAGS MATCHES
     "${QRAFT_NATIVE_CPU_FLAG_PATTERN}")
    message(FATAL_ERROR
      "QRAFT_STATIC_EXECUTABLE found native CPU tuning in configured "
      "compiler or linker flags")
  endif()
  if(NOT QRAFT_USE_JEMALLOC AND NOT WIN32)
    message(FATAL_ERROR
      "QRAFT_STATIC_EXECUTABLE requires static jemalloc. Install the "
      "jemalloc development package, or set QRAFT_STATIC_EXECUTABLE=OFF "
      "for a non-production developer build.")
  endif()
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(TF_BUILD_TESTS OFF CACHE BOOL "Build Taskflow tests" FORCE)
set(TF_BUILD_EXAMPLES OFF CACHE BOOL "Build Taskflow examples" FORCE)

add_subdirectory(external/taskflow)

include(cmake/PrepareDeterministicTaskflow.cmake)
qraft_prepare_deterministic_taskflow(
  "${CMAKE_CURRENT_SOURCE_DIR}/external/taskflow"
  "${CMAKE_CURRENT_BINARY_DIR}/generated/taskflow-deterministic")
target_include_directories(Taskflow BEFORE INTERFACE
  "$<BUILD_INTERFACE:${QRAFT_TASKFLOW_INCLUDE_ROOT}>")

find_package(Threads REQUIRED)

set(QRAFT_SOURCES
  src/core/bitset.cpp
  src/core/newick.cpp
  src/core/taxa.cpp
  src/core/tree.cpp
  src/genes/gene_index.cpp
  src/genes/gene_lca.cpp
  src/mre/mre.cpp
  src/mre/exact_split_ids.cpp
  src/mre/split_records.cpp
  src/mre/starting_tree.cpp
  src/mre/guide_tree.cpp
  src/counting/counters.cpp
  src/counting/units.cpp
  src/counting/hdt.cpp
  src/counting/extract_contract.cpp
  src/counting/count.cpp
  src/counting/query_taxa.cpp
  src/counting/query_tree.cpp
  src/divide/divide.cpp
  src/divide/bipartition_scores.cpp
  src/divide/bipartition_test.cpp
  src/regraft/regraft.cpp
  src/nni/nni.cpp
  src/nni/skeleton.cpp
  src/nni/spr.cpp
  src/solver/infer.cpp
  src/solver/support.cpp
  src/solver/recursion.cpp
  src/solver/species_tree.cpp
  src/solver/stats.cpp
  src/util/allocator.cpp
  src/util/parallel.cpp
  src/util/progress.cpp
)

add_library(qraft_lib ${QRAFT_SOURCES})
set_target_properties(qraft_lib PROPERTIES OUTPUT_NAME qraft)

target_include_directories(qraft_lib
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_compile_options(qraft_lib PRIVATE
  $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang>:-Wall>
  $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang>:-Wextra>
  $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang>:-Wpedantic>
)

if(QRAFT_NATIVE_OPTIMIZATIONS)
  target_compile_options(qraft_lib PRIVATE
    $<$<AND:$<CONFIG:Release>,$<COMPILE_LANG_AND_ID:CXX,GNU,Clang>>:-march=native>
  )
endif()

target_link_libraries(qraft_lib PUBLIC Threads::Threads PRIVATE Taskflow)

set(QRAFT_ALL_COMPILER_FLAGS
    "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS_RELEASE}")
if(WIN32)
  message(STATUS "Windows: using the system allocator")
elseif(QRAFT_STATIC_EXECUTABLE)
  if(QRAFT_ALL_COMPILER_FLAGS MATCHES "sanitize")
    message(FATAL_ERROR
      "Sanitizers are incompatible with the default full-static jemalloc "
      "profile. Configure the sanitizer build with "
      "QRAFT_STATIC_EXECUTABLE=OFF.")
  endif()

  set(QRAFT_SAVED_FIND_LIBRARY_SUFFIXES "${CMAKE_FIND_LIBRARY_SUFFIXES}")
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  find_library(QRAFT_JEMALLOC_STATIC_ARCHIVE NAMES jemalloc)
  set(CMAKE_FIND_LIBRARY_SUFFIXES "${QRAFT_SAVED_FIND_LIBRARY_SUFFIXES}")
  unset(QRAFT_SAVED_FIND_LIBRARY_SUFFIXES)
  if(NOT QRAFT_JEMALLOC_STATIC_ARCHIVE OR
     NOT EXISTS "${QRAFT_JEMALLOC_STATIC_ARCHIVE}" OR
     NOT QRAFT_JEMALLOC_STATIC_ARCHIVE MATCHES "\\.a$")
    message(FATAL_ERROR
      "The default static build requires libjemalloc.a, but CMake could not "
      "find it. Install libjemalloc-dev (Debian/Ubuntu) or the equivalent "
      "jemalloc development package.")
  endif()

  target_link_libraries(qraft_lib PUBLIC "${QRAFT_JEMALLOC_STATIC_ARCHIVE}")
  target_compile_definitions(qraft_lib PRIVATE QRAFT_LINKED_JEMALLOC=1)
  message(STATUS "Using static jemalloc: ${QRAFT_JEMALLOC_STATIC_ARCHIVE}")
elseif(QRAFT_USE_JEMALLOC AND NOT QRAFT_ALL_COMPILER_FLAGS MATCHES "sanitize")
  find_library(QRAFT_JEMALLOC_LIBRARY NAMES jemalloc jemalloc.so.2)
  if(NOT QRAFT_JEMALLOC_LIBRARY)
    find_file(QRAFT_JEMALLOC_LIBRARY NAMES libjemalloc.so.2
      PATHS /lib /usr/lib /usr/local/lib
      PATH_SUFFIXES x86_64-linux-gnu aarch64-linux-gnu)
  endif()
  if(QRAFT_JEMALLOC_LIBRARY)
    target_link_libraries(qraft_lib PUBLIC ${QRAFT_JEMALLOC_LIBRARY})
    target_compile_definitions(qraft_lib PRIVATE QRAFT_LINKED_JEMALLOC=1)
    message(STATUS "Using jemalloc: ${QRAFT_JEMALLOC_LIBRARY}")
  else()
    message(STATUS "jemalloc not found; using the platform allocator")
  endif()
endif()

add_executable(qraft
  src/cli/main.cpp
  src/cli/cli_options.cpp
  src/cli/gene_input.cpp)
target_link_libraries(qraft PRIVATE qraft_lib)

# qraft --version: the version, the git commit and the build options, written
# at build time so the commit stays current (cmake/WriteVersion.cmake).
if(QRAFT_STATIC_EXECUTABLE)
  set(QRAFT_BUILD_LINKAGE "static")
else()
  set(QRAFT_BUILD_LINKAGE "dynamic")
endif()
if(QRAFT_JEMALLOC_STATIC_ARCHIVE OR QRAFT_JEMALLOC_LIBRARY)
  set(QRAFT_BUILD_ALLOCATOR "jemalloc")
else()
  set(QRAFT_BUILD_ALLOCATOR "system allocator")
endif()
if(QRAFT_NATIVE_OPTIMIZATIONS)
  set(QRAFT_BUILD_CPU "tuned for the build machine (-march=native)")
else()
  set(QRAFT_BUILD_CPU "generic ${CMAKE_SYSTEM_PROCESSOR}")
endif()
set(QRAFT_BUILD_INFO
  "${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}, ${CMAKE_BUILD_TYPE}, ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR}, ${QRAFT_BUILD_LINKAGE}, ${QRAFT_BUILD_ALLOCATOR}, ${QRAFT_BUILD_CPU}")
set(QRAFT_VERSION_HEADER "${CMAKE_CURRENT_BINARY_DIR}/generated/version/qraft_version.hpp")
add_custom_target(qraft_version_header
  COMMAND ${CMAKE_COMMAND}
          -DQRAFT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
          -DQRAFT_OUTPUT=${QRAFT_VERSION_HEADER}
          -DQRAFT_VERSION=${PROJECT_VERSION}
          "-DQRAFT_BUILD_INFO=${QRAFT_BUILD_INFO}"
          -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/WriteVersion.cmake
  BYPRODUCTS ${QRAFT_VERSION_HEADER}
  COMMENT "Stamping the qraft version"
  VERBATIM)
add_dependencies(qraft qraft_version_header)
target_include_directories(qraft PRIVATE
  "${CMAKE_CURRENT_BINARY_DIR}/generated/version")
if(QRAFT_STATIC_EXECUTABLE)
  include(CheckLinkerFlag)
  check_linker_flag(CXX "-static" QRAFT_LINKER_SUPPORTS_FULL_STATIC)
  if(NOT QRAFT_LINKER_SUPPORTS_FULL_STATIC)
    message(FATAL_ERROR
      "The selected C++ toolchain cannot link a fully static executable")
  endif()
  target_link_options(qraft PRIVATE -static)
  if(WIN32)
    # No link time in the PE header, so a rebuild gives the same file.
    target_link_options(qraft PRIVATE -Wl,--no-insert-timestamp)
  endif()
endif()

if(QRAFT_ENABLE_IPO AND CMAKE_BUILD_TYPE STREQUAL "Release")
  include(CheckIPOSupported)
  check_ipo_supported(RESULT QRAFT_IPO_SUPPORTED OUTPUT QRAFT_IPO_ERROR)
  if(QRAFT_IPO_SUPPORTED)
    set_property(TARGET qraft_lib PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    set_property(TARGET qraft PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  else()
    message(STATUS "IPO/LTO not enabled: ${QRAFT_IPO_ERROR}")
  endif()
endif()

if(QRAFT_BUILD_TESTS)
  enable_testing()

  if(NOT QRAFT_NATIVE_OPTIMIZATIONS)
    string(CONCAT QRAFT_COMPILE_FLAG_MANIFEST_CONTENT
      "library=$<JOIN:$<TARGET_PROPERTY:qraft_lib,COMPILE_OPTIONS>, >\n"
      "library_definitions="
      "$<JOIN:$<TARGET_PROPERTY:qraft_lib,COMPILE_DEFINITIONS>, >\n"
      "cli=$<JOIN:$<TARGET_PROPERTY:qraft,COMPILE_OPTIONS>, >\n"
      "link=$<JOIN:$<TARGET_PROPERTY:qraft,LINK_OPTIONS>, >\n"
      "configured=${QRAFT_CONFIGURED_CXX_AND_LINK_FLAGS}\n"
      "debug=$<$<CONFIG:Debug>:${CMAKE_CXX_FLAGS_DEBUG}>\n"
      "release=$<$<CONFIG:Release>:${CMAKE_CXX_FLAGS_RELEASE}>\n"
      "relwithdebinfo="
      "$<$<CONFIG:RelWithDebInfo>:${CMAKE_CXX_FLAGS_RELWITHDEBINFO}>\n"
      "minsizerel=$<$<CONFIG:MinSizeRel>:${CMAKE_CXX_FLAGS_MINSIZEREL}>\n")
    file(GENERATE
      OUTPUT
        "${CMAKE_CURRENT_BINARY_DIR}/qraft-compile-flags-$<CONFIG>.txt"
      CONTENT "${QRAFT_COMPILE_FLAG_MANIFEST_CONTENT}"
      TARGET qraft_lib)
    add_test(
      NAME qraft.generic_isa_contract
      COMMAND ${CMAKE_COMMAND}
              -DQRAFT_FLAG_MANIFEST=${CMAKE_CURRENT_BINARY_DIR}/qraft-compile-flags-$<CONFIG>.txt
              -DQRAFT_LIBRARY=$<TARGET_FILE:qraft_lib>
              -DQRAFT_NM=${CMAKE_NM}
              -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/generic_isa_contract.cmake)
    set_tests_properties(qraft.generic_isa_contract PROPERTIES
      LABELS "regression;portability" TIMEOUT 60)
  endif()

  if(QRAFT_STATIC_EXECUTABLE AND NOT WIN32)
    set(QRAFT_READELF_EXECUTABLE "${CMAKE_READELF}")
    if(NOT QRAFT_READELF_EXECUTABLE OR
       NOT EXISTS "${QRAFT_READELF_EXECUTABLE}")
      find_program(QRAFT_READELF_EXECUTABLE NAMES readelf llvm-readelf)
    endif()
    if(NOT QRAFT_READELF_EXECUTABLE)
      message(FATAL_ERROR
        "readelf or llvm-readelf is required to test the static executable")
    endif()
    add_test(
      NAME qraft.static_cli_linkage
      COMMAND ${CMAKE_COMMAND}
              -DQRAFT_CLI=$<TARGET_FILE:qraft>
              -DQRAFT_READELF=${QRAFT_READELF_EXECUTABLE}
              -DQRAFT_NM=${CMAKE_NM}
              -DQRAFT_EXPECT_JEMALLOC=ON
              -DQRAFT_CROSSCOMPILING=${CMAKE_CROSSCOMPILING}
              -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/static_cli_linkage.cmake)
    set_tests_properties(qraft.static_cli_linkage PROPERTIES
      LABELS "regression;portability;packaging" TIMEOUT 60)
  endif()

  add_library(qraft_test_support STATIC tests/support/rf.cpp)
  target_include_directories(qraft_test_support
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/tests)
  target_link_libraries(qraft_test_support PUBLIC qraft_lib)

  add_executable(qraft_tests tests/qraft_tests.cpp)
  target_include_directories(qraft_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
  target_compile_definitions(qraft_tests PRIVATE
    QRAFT_TEST_TREE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/trees")
  target_link_libraries(qraft_tests PRIVATE qraft_test_support qraft_lib Taskflow)

  add_executable(qraft_recursion_tests
    tests/recursion_tests.cpp)
  target_include_directories(qraft_recursion_tests
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
  target_link_libraries(qraft_recursion_tests
    PRIVATE qraft_test_support qraft_lib Taskflow)
  foreach(test_case IN ITEMS
      parallel_loop_exception_safety
      parallel_child_exception_determinism
      concurrent_inference
      gene_restriction_duplicates
      recursion_complexity
      deep_skewed_recursion
      recursion_failure_unwinding)
    add_test(NAME qraft.${test_case}
      COMMAND qraft_recursion_tests ${test_case})
  endforeach()
  set_tests_properties(
    qraft.parallel_loop_exception_safety
    qraft.parallel_child_exception_determinism
    PROPERTIES LABELS "regression;safety" TIMEOUT 60)
  set_tests_properties(qraft.concurrent_inference PROPERTIES
    LABELS "regression;safety;determinism" TIMEOUT 60)
  set_tests_properties(qraft.gene_restriction_duplicates
    PROPERTIES LABELS "regression;compression;semantic" TIMEOUT 60)
  set_tests_properties(qraft.recursion_complexity PROPERTIES
    LABELS "regression;complexity" TIMEOUT 60)
  set_tests_properties(qraft.deep_skewed_recursion PROPERTIES
    LABELS "regression;safety;complexity" TIMEOUT 300)
  set_tests_properties(qraft.recursion_failure_unwinding PROPERTIES
    LABELS "regression;safety;determinism" TIMEOUT 300)

  add_test(
    NAME qraft.cli_parallel_parse_failure
    COMMAND ${CMAKE_COMMAND}
            -DQRAFT_CLI=$<TARGET_FILE:qraft>
            -DQRAFT_CLI_TEST_DIR=${CMAKE_CURRENT_BINARY_DIR}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/cli_parallel_parse_failure.cmake)
  set_tests_properties(qraft.cli_parallel_parse_failure PROPERTIES
    LABELS "regression;safety" TIMEOUT 60)

  add_test(
    NAME qraft.cli_spr
    COMMAND ${CMAKE_COMMAND}
            -DQRAFT_CLI=$<TARGET_FILE:qraft>
            -DQRAFT_CLI_TEST_DIR=${CMAKE_CURRENT_BINARY_DIR}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/cli_spr.cmake)
  set_tests_properties(qraft.cli_spr PROPERTIES
    LABELS "regression;determinism;dispatch" TIMEOUT 120)
  add_test(
    NAME qraft.cli_front
    COMMAND ${CMAKE_COMMAND}
            -DQRAFT_CLI=$<TARGET_FILE:qraft>
            -DQRAFT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
            -DQRAFT_CLI_TEST_DIR=${CMAKE_CURRENT_BINARY_DIR}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/cli_front.cmake)
  set_tests_properties(qraft.cli_front PROPERTIES
    LABELS "regression;dispatch" TIMEOUT 600)

  add_test(
    NAME qraft.cli_ablations
    COMMAND ${CMAKE_COMMAND}
            -DQRAFT_CLI=$<TARGET_FILE:qraft>
            -DQRAFT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
            -DQRAFT_CLI_TEST_DIR=${CMAKE_CURRENT_BINARY_DIR}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/cli_ablations.cmake)
  set_tests_properties(qraft.cli_ablations PROPERTIES
    LABELS "regression;determinism;dispatch" TIMEOUT 2400)

  if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    add_test(
      NAME qraft.cli_deep_inputs
      COMMAND ${CMAKE_COMMAND}
              -DQRAFT_CLI=$<TARGET_FILE:qraft>
              -DQRAFT_CLI_TEST_DIR=${CMAKE_CURRENT_BINARY_DIR}
              -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/cli_deep_inputs.cmake)
    set_tests_properties(qraft.cli_deep_inputs PROPERTIES
      LABELS "regression;safety;determinism" TIMEOUT 600)
  endif()

  add_test(
    NAME qraft.deterministic_containers
    COMMAND ${CMAKE_COMMAND}
            -DQRAFT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
            -DQRAFT_TASKFLOW_EXECUTOR=${QRAFT_TASKFLOW_EXECUTOR}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_deterministic_containers.cmake)
  set_tests_properties(qraft.deterministic_containers PROPERTIES
    LABELS "proof;complexity;determinism" TIMEOUT 60)

  set(QRAFT_BASIC_TEST_CASES
    newick_parser
    newick_structure
    newick_deep_trees
    gene_structure_validation
    checked_arithmetic_guards
    hdt_fail_closed_guards
    taxon_registry_scaling
    mre_duplicate_and_incomplete_genes
    mre_exact_cluster_ids
    mre_exact_cluster_ids_parallel
    mre_exact_split_entries
    mre_update_stream_stop
    incomplete_gene_trees
    trivial_bipartition_completion
    bipartition_test
    bipartition_test_distance_index
    bipartition_test_incomplete_genes
    bipartition_test_shared_genes_reduce_to_j
    bipartition_test_shared_genes_noise_free
    nni_split_support
    nni_exact_split_support
    bipartition_scores_match_enumeration
    regraft_scores_match_enumeration
    skeleton_gains_exactness
    spr_exactness
    spr_polytomy_exactness
    spr_gain_exactness
    guide_tree
    guide_seed
    nni_gains_match_enumeration
    branch_support_matches_enumeration
    branch_support_output
    nni_binary_rooted_starting_tree
    nni_local_pass_matches_full_pass
    nni_batch_gains_add
    wide_count_exact
    bipartition_test_threshold_values
    msc_simulation_recovers_species_tree
    child_gene_set_partition_bound
    mre_complexity_primitives
    counters_packed_color_masks
    counters_match_enumeration
    hdt_units
    hdt_construction
    hdt_contract
    hdt_extract
    count_traversal_oracle
    count_traversal_scaling
    qraft_smoke
    worker_cap_by_input_size
    worker_count_determinism
    pipeline_determinism
    ablation_smoke
    hdt_weighted_fixtures
    hdt_duplicate_owner_contexts
    tie_band_selection_permutations
    bipartition_near_ties
    regraft_near_ties
    hdt_heavy_threshold
    hdt_structural_bounds
    finite_score_guards
    rf_identity
    fixture_smoke
  )
  foreach(test_case IN LISTS QRAFT_BASIC_TEST_CASES)
    add_test(NAME qraft.${test_case} COMMAND qraft_tests ${test_case})
  endforeach()
  set_tests_properties(
    qraft.newick_structure
    qraft.newick_deep_trees
    qraft.mre_update_stream_stop
    PROPERTIES LABELS "regression;safety" TIMEOUT 300)

  find_package(Python3 COMPONENTS Interpreter QUIET)
  if(Python3_Interpreter_FOUND)
    # Exact state/transition closure, quotient factorization,
    # parenthesization, and corridor-orientation laws.
    add_test(
      NAME qraft.counting_scheme_exhaustive
      COMMAND ${Python3_EXECUTABLE} -B
              ${CMAKE_CURRENT_SOURCE_DIR}/tools/counting/check_counting_scheme.py
    )
    add_test(
      NAME qraft.counters_generated
      COMMAND ${Python3_EXECUTABLE} -B
              ${CMAKE_CURRENT_SOURCE_DIR}/tools/counting/generate_counters.py
              --check-full
    )
    add_test(
      NAME qraft.counters_reproducible
      COMMAND ${Python3_EXECUTABLE} -B
              ${CMAKE_CURRENT_SOURCE_DIR}/tools/counting/generate_counters.py
              --check-reproducible
    )

    # The machine checks of the counting scheme (verification/) and the check that
    # the production tables are its natural counting scheme.
    foreach(check IN ITEMS production_counters natural_closure span_full
                           end_to_end lower_bound)
      add_test(
        NAME qraft.verification_${check}
        COMMAND ${Python3_EXECUTABLE} -B
                ${CMAKE_CURRENT_SOURCE_DIR}/verification/check_${check}.py
      )
      set_tests_properties(qraft.verification_${check} PROPERTIES
        LABELS "proof;verification" TIMEOUT 300)
    endforeach()

    # Deterministic all-topology/all-contraction reference shard. The full
    # through-eight command uses the same seed and 256 shards.
    add_test(
      NAME qraft.small_trees_shard_0_of_8
      COMMAND ${Python3_EXECUTABLE} -B
              ${CMAKE_CURRENT_SOURCE_DIR}/tools/counting/check_small_trees.py
              --max-leaves 6 --shard-count 8 --shard-index 0
              --seed 20260713
    )

    # The actual production quotient/base/merge/Contract/ExtractContract path.
    # Two small deterministic shards stay in the normal suite.  Configure with
    # QRAFT_REGISTER_LONG_TESTS=ON to register the complete through-eight set.
    set(QRAFT_PRODUCTION_SMALL_REPRESENTATIVE_TESTS)
    foreach(production_shard IN ITEMS 0 37)
      set(production_test
        qraft.production_small_trees_shard_${production_shard}_of_64)
      add_test(
        NAME ${production_test}
        COMMAND ${Python3_EXECUTABLE} -B
                ${CMAKE_CURRENT_SOURCE_DIR}/tools/counting/check_production_small_trees.py
                --runner $<TARGET_FILE:qraft_tests>
                --max-leaves 6 --shard-count 64
                --shard-index ${production_shard} --seed 20260713
      )
      list(APPEND QRAFT_PRODUCTION_SMALL_REPRESENTATIVE_TESTS
        ${production_test})
    endforeach()

    if(QRAFT_REGISTER_LONG_TESTS)
      set(QRAFT_PRODUCTION_SMALL_LONG_TESTS)
      foreach(production_shard RANGE 0 255)
        set(production_test
          qraft.production_small_trees_through_8_shard_${production_shard}_of_256)
        add_test(
          NAME ${production_test}
          COMMAND ${Python3_EXECUTABLE} -B
                  ${CMAKE_CURRENT_SOURCE_DIR}/tools/counting/check_production_small_trees.py
                  --runner $<TARGET_FILE:qraft_tests>
                  --max-leaves 8 --shard-count 256
                  --shard-index ${production_shard} --seed 20260713
        )
        list(APPEND QRAFT_PRODUCTION_SMALL_LONG_TESTS ${production_test})
      endforeach()
      set_tests_properties(${QRAFT_PRODUCTION_SMALL_LONG_TESTS} PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1;PYTHONHASHSEED=0"
        LABELS "proof;semantic;exhaustive;long"
        TIMEOUT 1800)
    endif()

    # Independent quartet-count negative controls.
    add_test(
      NAME qraft.quartet_oracle
      COMMAND ${Python3_EXECUTABLE} -B
              ${CMAKE_CURRENT_SOURCE_DIR}/tools/counting/check_quartets.py
    )

    set(QRAFT_PYTHON_PROOF_TESTS
      qraft.counting_scheme_exhaustive
      qraft.counters_generated
      qraft.counters_reproducible
      qraft.small_trees_shard_0_of_8
      ${QRAFT_PRODUCTION_SMALL_REPRESENTATIVE_TESTS}
      qraft.quartet_oracle
    )
    set_tests_properties(${QRAFT_PYTHON_PROOF_TESTS} PROPERTIES
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      ENVIRONMENT PYTHONDONTWRITEBYTECODE=1)
    set_tests_properties(
      qraft.counting_scheme_exhaustive
      qraft.small_trees_shard_0_of_8
      PROPERTIES LABELS "proof;exhaustive" TIMEOUT 300)
    set_tests_properties(
      ${QRAFT_PRODUCTION_SMALL_REPRESENTATIVE_TESTS}
      PROPERTIES
        ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1;PYTHONHASHSEED=0"
        LABELS "proof;semantic;exhaustive"
        TIMEOUT 300)
    set_tests_properties(
      qraft.counters_generated
      qraft.counters_reproducible
      PROPERTIES LABELS "proof;reproducibility" TIMEOUT 180)
    set_tests_properties(qraft.quartet_oracle
      PROPERTIES LABELS "proof;semantic" TIMEOUT 60)
  endif()

  # Supplemental fixed-fixture production coverage for quotient DP, base HDT,
  # Contract, ExtractContract, and connected traversal contexts.
  set_tests_properties(
    qraft.counters_packed_color_masks
    qraft.counters_match_enumeration
    qraft.hdt_construction
    qraft.hdt_contract
    qraft.hdt_extract
    qraft.count_traversal_oracle
    PROPERTIES LABELS "proof;semantic")
  set_tests_properties(
    qraft.finite_score_guards
    qraft.gene_structure_validation
    qraft.checked_arithmetic_guards
    qraft.hdt_fail_closed_guards
    PROPERTIES LABELS "proof;assumption;safety")

  # Named semantic and structural obligations from Complexity.md.
  set_tests_properties(
    qraft.hdt_weighted_fixtures
    qraft.hdt_duplicate_owner_contexts
    qraft.hdt_heavy_threshold
    qraft.bipartition_scores_match_enumeration
    qraft.regraft_scores_match_enumeration
    qraft.skeleton_gains_exactness
    qraft.spr_exactness
    qraft.spr_polytomy_exactness
    qraft.spr_gain_exactness
    qraft.nni_gains_match_enumeration
    qraft.branch_support_matches_enumeration
    qraft.nni_binary_rooted_starting_tree
    qraft.nni_local_pass_matches_full_pass
    qraft.nni_batch_gains_add
    qraft.wide_count_exact
    qraft.bipartition_test_threshold_values
    qraft.bipartition_test_shared_genes_reduce_to_j
    qraft.bipartition_test_shared_genes_noise_free
    qraft.msc_simulation_recovers_species_tree
    PROPERTIES LABELS "proof;semantic")
  set_tests_properties(
    qraft.taxon_registry_scaling
    qraft.mre_duplicate_and_incomplete_genes
    qraft.incomplete_gene_trees
    qraft.child_gene_set_partition_bound
    qraft.mre_complexity_primitives
    qraft.hdt_structural_bounds
    PROPERTIES LABELS "proof;structural")
  set_property(TEST qraft.taxon_registry_scaling APPEND
    PROPERTY LABELS "complexity;determinism")
  set_property(TEST qraft.mre_duplicate_and_incomplete_genes APPEND
    PROPERTY LABELS "determinism")
  set_property(TEST qraft.child_gene_set_partition_bound APPEND
    PROPERTY LABELS "complexity")
  set_property(TEST qraft.mre_complexity_primitives APPEND
    PROPERTY LABELS "complexity;determinism")
  set_property(TEST qraft.worker_cap_by_input_size APPEND
    PROPERTY LABELS "complexity;determinism;safety")

  # Powers through 32768 for every documented structural family. It is labeled
  # large so local fast runs can use `ctest -LE large`.
  add_test(NAME qraft.hdt_structural_bounds_32768
           COMMAND qraft_tests hdt_structural_bounds)
  set_tests_properties(qraft.hdt_structural_bounds_32768 PROPERTIES
    ENVIRONMENT QRAFT_TEST_STRUCTURAL_MAX=32768
    LABELS "proof;structural;large"
    TIMEOUT 300)

  # Operation-count and measured high-water gates.
  set_tests_properties(qraft.count_traversal_scaling PROPERTIES
    LABELS "proof;complexity;memory")
  set_property(TEST qraft.hdt_heavy_threshold APPEND PROPERTY
    LABELS complexity)

  # Clean worker-count determinism for active split, one-X, and NNI public
  # workflows.
  set_tests_properties(
    qraft.worker_count_determinism
    qraft.pipeline_determinism
    qraft.trivial_bipartition_completion
    qraft.bipartition_test
    qraft.nni_split_support
    qraft.tie_band_selection_permutations
    qraft.bipartition_near_ties
    qraft.regraft_near_ties
    PROPERTIES LABELS "proof;determinism")
endif()
