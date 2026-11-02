# QRAFT

QRAFT (Quartet Regrafting) infers a species tree from gene trees. It looks for
the species tree that agrees with the most gene-tree quartets, counting every
quartet exactly, and it scales to very large datasets. The tree it returns does
not depend on the number of threads.

- Input: gene trees, one Newick tree per line. Gene trees may miss taxa or
  hold several copies of a taxon.
- Output: the species tree as one Newick line, optionally with local posterior
  probabilities on its edges.

## Install

Download the release archive for Linux x86-64, unpack it, and check it:

```bash
tar xf qraft-1.0.0-linux-x86_64.tar.gz
cd qraft-1.0.0-linux-x86_64
./install-check.sh          # runs the built-in selftest and the example
export PATH="$PWD/bin:$PATH"
```

The binary is static: it needs no libraries and runs on any x86-64 Linux, and
under WSL2 on Windows.

To build from source (CMake 3.20 or newer, a C++20 compiler, and jemalloc's
static library, `libjemalloc-dev` on Debian and Ubuntu):

```bash
git clone --recursive https://github.com/risenfromashes/qraft.git
cd qraft
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build -R '^qraft\.' -j 8
```

Without jemalloc, configure with `-DQRAFT_STATIC_EXECUTABLE=OFF -DQRAFT_USE_JEMALLOC=OFF`.

## Quick start

```bash
qraft -i genes.tre -o species.tre                 # infer a species tree
qraft -i genes.tre -o species.tre -t 8            # on 8 threads
qraft -i genes.tre -o species.tre --mode extended # slower, for noisy data
qraft -i genes.tre -o species.tre --support localpp
qraft support -i genes.tre --tree species.tre -o annotated.tre
qraft check genes.tre                             # summarize the input first
```

`examples/` holds a small simulated dataset (100 taxa, 200 gene trees) and the
tree this release infers from it:

```bash
qraft -i examples/genes-100taxa.tre -o out.tre
cmp out.tre examples/expected-100taxa.tre && echo same
```

## Commands

| command | what it does |
| --- | --- |
| `qraft -i GENES` (or `qraft infer`) | infer a species tree |
| `qraft support -i GENES --tree TREE` | write local posterior probabilities on a given tree |
| `qraft check GENES` | summarize a gene-tree file: taxa, missing and repeated taxa, polytomies |
| `qraft version` | version, commit and build details |
| `qraft selftest` | run a small built-in example and check the result |
| `qraft completion bash\|zsh\|fish` | print a shell completion script |
| `qraft help TOPIC` | options, modes, guide, support, ablations, test, output |

Every option is listed by `qraft help options` and in the man page
(`man/man1/qraft.1` in the release archive).

## Modes and options

- **Default mode**: the majority-rule extended consensus of the gene trees,
  a recursive divide along well-supported bipartitions with a statistical test
  of each bipartition, contextual regrafting of the solved sides, and
  nearest-neighbour interchanges.
- **Extended mode** (`--mode extended`): the same, then unbounded batched
  NNI and an exact SPR refinement. It is slower and meant for noisy biological
  data.
- **Guide tree** (`--guide-tree PATH`): edits the starting tree with a trusted
  tree on some or all taxa before the run.
- **Branch support** (`--support localpp`): local posterior probability and
  length in coalescent units on every internal edge whose ends are binary.

Progress and a run summary go to standard error on a terminal; `-q` prints
errors only and `-v` prints the summary in scripts too. `--stats-json PATH`
writes the options, the version and run statistics. Exit status: 0 on
success, 1 when an input cannot be used or the run fails, 2 for a command-line
error, 3 when memory runs out.

To compare the quartet scores of given trees, use tqDist (Sand et al.,
Bioinformatics, 2014).

## Citing QRAFT

Ashrafur Rahman, Anik Saha and Md. Shamsuzzoha Bayzid. QRAFT: Scaling
coalescent-based species tree inference to one million taxa. Manuscript in
preparation.

`CITATION.cff` holds the same reference in machine-readable form.

## License

Copyright (C) 2026 Ashrafur Rahman, Anik Saha and Md. Shamsuzzoha Bayzid.

QRAFT is free software under the GNU Affero General Public License, version 3
or (at your option) any later version; see `LICENSE`. The release binary
includes third-party code under its own licenses; see
`THIRD_PARTY_NOTICES.md`.
