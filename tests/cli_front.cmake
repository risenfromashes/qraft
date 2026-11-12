# The command-line front: commands, help topics, error messages and exit
# codes, every option accepted, help/man/completion consistency, and the new
# command forms giving byte-identical trees to the long-flag forms.

cmake_policy(VERSION 3.20)

if(NOT DEFINED QRAFT_CLI OR NOT DEFINED QRAFT_SOURCE_DIR OR
   NOT DEFINED QRAFT_CLI_TEST_DIR)
  message(FATAL_ERROR "QRAFT_CLI, QRAFT_SOURCE_DIR, and QRAFT_CLI_TEST_DIR are required")
endif()

set(dir "${QRAFT_CLI_TEST_DIR}/cli-front")
file(REMOVE_RECURSE "${dir}")
file(MAKE_DIRECTORY "${dir}")
set(genes "${QRAFT_SOURCE_DIR}/tests/trees/astral2_200taxa_50genes.tre")
set(missing_genes "${QRAFT_SOURCE_DIR}/tests/trees/astral2_200taxa_50genes_25missing.tre")

# Runs the CLI; sets out, err and code in the caller.
function(run)
  execute_process(COMMAND "${QRAFT_CLI}" ${ARGN}
                  WORKING_DIRECTORY "${dir}"
                  RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
  set(out "${stdout}" PARENT_SCOPE)
  set(err "${stderr}" PARENT_SCOPE)
  set(code "${result}" PARENT_SCOPE)
endfunction()

# Expects exit status `status`, and `pattern` in standard output (stream out)
# or standard error (stream err).
function(expect status stream pattern)
  run(${ARGN})
  if(stream STREQUAL "out")
    set(text "${out}")
  else()
    set(text "${err}")
  endif()
  if(NOT code EQUAL status OR NOT text MATCHES "${pattern}")
    message(FATAL_ERROR "qraft ${ARGN}: expected status ${status} and /${pattern}/ "
                        "on ${stream}; got ${code}\nstdout:\n${out}\nstderr:\n${err}")
  endif()
endfunction()

function(same_file a b what)
  file(SHA256 "${a}" ha)
  file(SHA256 "${b}" hb)
  if(NOT ha STREQUAL hb)
    message(FATAL_ERROR "${what}: ${a} and ${b} differ")
  endif()
endfunction()

# Version and help.
expect(0 out "^qraft [0-9]+\\.[0-9]+\\.[0-9]+ \\(([0-9a-f]+(-dirty)?|unknown)\\)" --version)
expect(0 out "license: AGPL-3.0-or-later" version)
expect(0 out "build: " version)
expect(0 out "Commands:.*infer.*support.*check.*Examples:" --help)
expect(0 out "Commands:" help)
foreach(topic options modes guide support ablations test output infer support check
              version selftest completion help)
  expect(0 out "." help ${topic})
endforeach()
expect(0 out "Usage: qraft support" support --help)
expect(0 out "Usage: qraft check" check --help)
expect(2 err "unknown help topic: optons \\(did you mean options\\?\\)" help optons)

# Help, man page and completions list the same options.
run(help options)
set(help_text "${out}")
string(REGEX MATCHALL "--[a-z][a-z0-9-]*" help_options "${help_text}")
list(REMOVE_DUPLICATES help_options)
run(help man)
string(REPLACE "\\-" "-" man_text "${out}")
run(completion bash)
set(bash_text "${out}")
run(completion fish)
set(fish_text "${out}")
foreach(option IN LISTS help_options)
  string(REPLACE "-" "\\-" escaped "${option}")
  if(NOT man_text MATCHES "${escaped}")
    message(FATAL_ERROR "the man page lacks ${option}")
  endif()
  if(NOT bash_text MATCHES "${escaped}")
    message(FATAL_ERROR "the bash completion lacks ${option}")
  endif()
  string(SUBSTRING "${option}" 2 -1 bare)
  if(NOT fish_text MATCHES "-l ${bare}( |\n)")
    message(FATAL_ERROR "the fish completion lacks ${option}")
  endif()
endforeach()
expect(0 out "bashcompinit" completion zsh)
expect(2 err "unknown shell: tcsh" completion tcsh)
# The completion scripts parse in the shells that are installed.
foreach(shell bash zsh fish)
  find_program(shell_program_${shell} ${shell})
  set(shell_program "${shell_program_${shell}}")
  if(NOT shell_program)
    continue()
  endif()
  run(completion ${shell})
  file(WRITE "${dir}/completion.${shell}" "${out}")
  execute_process(COMMAND "${shell_program}" -n "${dir}/completion.${shell}"
                  RESULT_VARIABLE syntax ERROR_VARIABLE syntax_err)
  if(NOT syntax EQUAL 0)
    message(FATAL_ERROR "the ${shell} completion does not parse: ${syntax_err}")
  endif()
endforeach()

# Selftest.
expect(0 out "selftest passed" selftest -t 2)
expect(2 err "invalid value for -t: x \\(expected a whole number\\)" selftest -t x)
expect(2 err "missing value for -t" selftest -t)

# Errors, with plain messages and exit codes.
expect(1 err "cannot read the gene trees file 'nosuch.tre': no such file" -i nosuch.tre)
file(WRITE "${dir}/empty.tre" "")
expect(1 err "the gene trees file 'empty.tre' is empty" -i empty.tre)
expect(1 err "the gene trees file '/dev/null' is a pipe or special file" -i /dev/null)
file(WRITE "${dir}/blank.tre" "\n  \n")
expect(1 err "holds no Newick trees" -i blank.tre)
file(WRITE "${dir}/bad.tre" "((A,B),(C,D));\n\n((A,B),(C,D);\n")
expect(1 err "line 3, is not a valid Newick tree: .*expected" -i bad.tre)
expect(2 err "unknown argument: --no-such-option" -i "${genes}" --no-such-option)
expect(2 err "unknown argument: --gene \\(did you mean --genes\\?\\)" --gene "${genes}")
expect(2 err "unknown command: chek \\(did you mean check\\?\\)" chek "${genes}")
expect(2 err "invalid value for --mode: fast \\(expected default\\|extended\\)"
       -i "${genes}" --mode fast)
expect(2 err "invalid value for --divide-score: best" -i "${genes}" --divide-score best)
expect(2 err "invalid value for --threads: x \\(expected a whole number\\)" -i "${genes}" -t x)
expect(2 err "missing value for --out" -i "${genes}" --out)
expect(2 err "--no-nni takes no value" -i "${genes}" --no-nni=1)
expect(2 err "--genes is required" -o out.tre)
expect(2 err "--tree is required" support -i "${genes}")
expect(2 err "--mode is not an option of 'qraft check'" check "${genes}" --mode extended)
expect(2 err "unexpected argument: extra.tre" -i "${genes}" extra.tre)
expect(2 err "--support-detail needs --support localpp" -i "${genes}" --support-detail)
expect(2 err "--guide-seed needs --guide-tree" -i "${genes}" --guide-seed)
file(WRITE "${dir}/some_guide.tre" "((A,B),(C,D));\n")
expect(2 err "--guide-tree and --starting-tree cannot be combined"
       -i "${genes}" --guide-tree some_guide.tre --starting-tree some_guide.tre)
expect(1 err "cannot write the output file '/nonexistent-qraft-dir/out.tre': the folder"
       -i "${genes}" -o /nonexistent-qraft-dir/out.tre)
expect(1 err "cannot read the guide tree file 'noguide.tre'" -i "${genes}" --guide-tree noguide.tre)
# Option values the library rejects are usage errors, found before reading.
expect(2 err "SPR refinement requires positive" -i "${genes}" --spr-rounds 0)
expect(2 err "bipartition test requires" -i "${genes}" --test-coefficient nan)
# Trees given as files are parsed before the gene trees are read.
file(WRITE "${dir}/bad_guide.tre" "((A,B),(C,D);\n")
expect(1 err "the guide tree file 'bad_guide.tre' is not a valid Newick tree"
       -i "${genes}" --guide-tree bad_guide.tre)
expect(1 err "the starting tree file 'bad_guide.tre' is not a valid Newick tree"
       -i "${genes}" --starting-tree bad_guide.tre)
# An output that is a link to a missing folder is found before the run.
execute_process(COMMAND ln -sf /nonexistent-qraft-dir/linked.tre "${dir}/dangling.tre")
expect(1 err "cannot write the output file 'dangling.tre': the folder '/nonexistent-qraft-dir'"
       -i "${genes}" -o dangling.tre)
expect(2 err "Usage:" )

# --guide-seed is a no-op kept for old command lines: a note, and the tree of
# plain --guide-tree.
run(-i "${genes}" -o noop_guide.tre -q -t 2)
run(-i "${genes}" -o noop_plain.tre -t 2 --guide-tree noop_guide.tre)
run(-i "${genes}" -o noop_seed.tre -t 2 --guide-tree noop_guide.tre --guide-seed)
if(NOT code EQUAL 0 OR NOT err MATCHES "--guide-seed has no effect since 1.0.1")
  message(FATAL_ERROR "--guide-seed: expected status 0 and the no-op note; got ${code}: ${err}")
endif()
same_file("${dir}/noop_plain.tre" "${dir}/noop_seed.tre" "--guide-seed is a no-op")

# Quiet runs print nothing on standard error, even with --progress.
run(-i "${genes}" -o quiet.tre -q -t 2)
if(NOT code EQUAL 0 OR NOT err STREQUAL "")
  message(FATAL_ERROR "-q printed on standard error: ${err}")
endif()
run(-i "${genes}" -o quiet_progress.tre -q --progress -t 2)
if(NOT code EQUAL 0 OR NOT err STREQUAL "")
  message(FATAL_ERROR "-q --progress printed on standard error: ${err}")
endif()
# Verbose runs print the start line and the summary.
expect(0 err "reading .* threads.*done in .* taxa, 50 gene trees" -i "${genes}" -o verbose.tre -v -t 2)

# check.
expect(0 out "50 gene trees on 201 taxa.*genes missing taxa: 0.*OK: " check "${genes}")
expect(0 out "genes missing taxa: [1-9].*note: genes miss taxa" check -i "${missing_genes}" --json check.json)
file(READ "${dir}/check.json" check_json)
string(JSON usable GET "${check_json}" usable)
string(JSON genes_count GET "${check_json}" genes)
if(NOT usable STREQUAL "ON" OR NOT genes_count EQUAL 50)
  message(FATAL_ERROR "check --json: ${check_json}")
endif()
file(WRITE "${dir}/tiny.tre" "(A,B,C);\n(A,B,C);\n")
expect(1 err "fewer than 4 taxa" check tiny.tre)

# The new forms give the trees of the long-flag forms.
run(--genes "${genes}" --out flags.tre --threads 2 --no-progress)
# A thread count below 1 means every hardware thread, as before.
run(-i "${genes}" -o all_threads.tre -t -1)
same_file("${dir}/flags.tre" "${dir}/all_threads.tre" "-t -1")
run(infer "${genes}" -o positional.tre -t 2)
same_file("${dir}/flags.tre" "${dir}/positional.tre" "qraft infer GENES")
run(--genes=${genes} --out=equals.tre --threads=3)
same_file("${dir}/flags.tre" "${dir}/equals.tre" "--option=value")
run(-i "${genes}" -o extended_flag.tre -t 2 --extended)
run(-i "${genes}" -o extended_mode.tre -t 2 --mode extended)
same_file("${dir}/extended_flag.tre" "${dir}/extended_mode.tre" "--mode extended")
run(-i "${genes}" -o reset.tre -t 2 --extended --mode default)
same_file("${dir}/flags.tre" "${dir}/reset.tre" "--mode default after --extended")
run(-i "${genes}" --starting-tree flags.tre --no-divide --no-nni --support localpp
    -o support_flags.tre -t 2)
run(support -i "${genes}" --tree flags.tre -o support_command.tre -t 2)
same_file("${dir}/support_flags.tre" "${dir}/support_command.tre" "qraft support")
run(support -i "${genes}" --tree flags.tre -o support_detail.tre -t 2 --support-detail
    --support-frequencies quartet)
if(NOT code EQUAL 0)
  message(FATAL_ERROR "qraft support --support-detail failed: ${err}")
endif()

# The statistics name the version.
run(-i "${genes}" -o stats.tre -t 2 --stats-json stats.json)
file(READ "${dir}/stats.json" stats)
string(JSON version GET "${stats}" qraft_version)
if(NOT version MATCHES "^[0-9]+\\.[0-9]+\\.[0-9]+ \\(")
  message(FATAL_ERROR "qraft_version missing from the statistics: ${stats}")
endif()

# Every option of the earlier command line is still accepted.
set(guide "${dir}/flags.tre")
foreach(case
    "--progress" "--no-progress" "--no-divide" "--divide-score;balanced"
    "--divide-score;quartet" "--no-bipartition-test" "--regraft-score;normalized"
    "--regraft-score;raw" "--no-nni" "--no-nni-cap" "--no-nni-support-gate"
    "--nni-batch" "--spr-refine" "--spr-reads;64" "--spr-min-clade;5"
    "--spr-candidates;20" "--spr-rounds;1" "--support;localpp" "--support;none"
    "--support;localpp;--support-frequencies;quartet"
    "--support;localpp;--support-frequencies;gene" "--support;localpp;--support-detail"
    "--test-coefficient;0.1" "--test-crossover;100" "--test-proposals;2"
    "--test-anchor-pairs;2" "--guide-tree;${guide}" "--guide-tree;${guide};--guide-seed"
    "--guide-tree;${guide};--write-starting-tree;written.tre"
    "--starting-tree;${guide}" "--stats-json;flags.json")
  run(-i "${genes}" -o accepted.tre -t 2 ${case})
  if(NOT code EQUAL 0)
    message(FATAL_ERROR "qraft ${case} failed with ${code}: ${err}")
  endif()
endforeach()
