# QRAFT


QRAFT (Quartet Regrafting) infers a species tree from gene trees. It looks for
the species tree that agrees with the most gene-tree quartets, counting every
quartet exactly, and it scales to very large datasets. The tree it returns does
not depend on the number of threads.

- Input: gene trees, one Newick tree per line. Gene trees may miss taxa or
  hold several copies of a taxon.
- Output: the species tree as one Newick line, optionally with local posterior
  probabilities on its edges.

On 64 CPU threads and no GPU, QRAFT analyzed one million taxa with 1,000 gene
trees in under 5.5 hours, and one million gene trees on 1,000 taxa in under 35
minutes.

## Install

Download the release archive for your system, unpack it, and check it.

Linux (x86-64):

```bash
tar xf qraft-1.0.0-linux-x86_64.tar.gz
cd qraft-1.0.0-linux-x86_64
./install-check.sh          # runs the built-in selftest and the example
export PATH="$PWD/bin:$PATH"
```

Windows (x86-64): unzip `qraft-1.0.0-windows-x86_64.zip`, then in a Command
Prompt in that folder:

```bat
install-check.cmd
qraft -i examples\genes-100taxa.tre -o species.tre
```

Both binaries are static: they need no libraries, and they infer the same trees.
The Linux binary uses jemalloc as its memory allocator and the Windows binary
the system allocator, so Windows runs may take longer. On Windows, keep file
paths to plain ASCII characters.

To build QRAFT yourself, see [Building from source](#building-from-source).

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

## Sample datasets

A sample of the large-scale dataset
([zip](https://huggingface.co/datasets/risenfromashes/qraft-data/resolve/main/datasets/scalability-sample/scalability-sample-10k-100k-1M-R1-R2.zip),
28 GB) holds two replicates each of the 10,000, 100,000 and 1,000,000-taxon
datasets, with 1,000 gene trees each. It can be downloaded and run as follows:

```bash
curl -LO https://huggingface.co/datasets/risenfromashes/qraft-data/resolve/main/datasets/scalability-sample/scalability-sample-10k-100k-1M-R1-R2.zip
unzip scalability-sample-10k-100k-1M-R1-R2.zip
qraft -i scalability-sample/t_10000_*/R1/all_gt.tre -o species-10k.tre
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
- **Refined starting tree** (`--refine-start`), the refined-MRE ablation of
  the paper: first refines the MRE starting tree as `--no-divide --extended`
  does, then runs from it. It gives the same trees as these runs:

```bash
qraft -i genes.tre -o mre.tre --no-divide --no-nni
qraft -i genes.tre -o refined.tre --starting-tree mre.tre --no-divide --extended
qraft -i genes.tre -o divided.tre --starting-tree refined.tre --no-nni  # = --refine-start --no-nni
qraft -i genes.tre -o out.tre --starting-tree divided.tre --no-divide   # = --refine-start
qraft -i genes.tre -o out.tre --starting-tree divided.tre --no-divide --extended  # = --refine-start --extended
```

`qraft help ablations` says which run each other option goes to.

Progress and a run summary go to standard error on a terminal; `-q` prints
errors only and `-v` prints the summary in scripts too. `--stats-json PATH`
writes the options, the version and run statistics. Exit status: 0 on
success, 1 when an input cannot be used or the run fails, 2 for a command-line
error, 3 when memory runs out.

To compare the quartet scores of given trees, use tqDist (Sand et al.,
Bioinformatics, 2014).

## Simulating the datasets

The paper's datasets were simulated with
[SimPhy-X](https://github.com/aaniksahaa/simphy-x) 1.0.0, a
multithreaded re-implementation of SimPhy 1.0.2 that gives the same trees for
the same command and seed. Two scripts in `simulation/` run it with the paper's
model and seed (Linux, bash, Python 3):

- `simulate.sh` simulates one dataset, given the number of taxa (`-t`), gene
  trees (`-g`) and replicates, and the population-size range (`--spmin`,
  `--spmax`) that sets the level of incomplete lineage sorting. Each replicate
  folder `R1`, `R2`, ... gets the gene trees in `all_gt.tre`, one per line, and
  the true species tree in `s_tree.trees`.
- `paper-datasets.sh` runs `simulate.sh` for every dataset of one experiment:
  `scalability`, `million` or `vary-ils` (`--dry-run` lists the commands).

```bash
curl -LO https://github.com/aaniksahaa/simphy-x/releases/download/v1.0.0/simphy-x-1.0.0-linux-x86_64
chmod +x simphy-x-1.0.0-linux-x86_64
# one dataset: 1,000 taxa, 1,000 gene trees, five replicates
simulation/simulate.sh -t 1000 -g 1000 --replicates 5 --spmin 100000 --spmax 200000 \
  --simphy-bin ./simphy-x-1.0.0-linux-x86_64
# all datasets of one experiment (up to 1.7 TB)
simulation/paper-datasets.sh scalability --simphy-bin ./simphy-x-1.0.0-linux-x86_64
```

## Building from source

QRAFT needs CMake 3.20 or newer and a C++20 compiler (tested with GCC 13 and
14). Its one dependency, Taskflow, comes as a git submodule.

Linux (x86-64), on Debian 13 or Ubuntu 24.04:

```bash
sudo apt-get update && sudo apt-get install -y g++ cmake git python3 libjemalloc-dev
git clone --recursive https://github.com/risenfromashes/qraft.git
cd qraft
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 4
ctest --test-dir build -R '^qraft\.' -j 4     # the test suite, a few minutes
build/qraft selftest
```

If you cloned without `--recursive`, fetch Taskflow with
`git submodule update --init external/taskflow`. The default build is a static
binary linked with jemalloc, so it needs jemalloc's static library
(`libjemalloc-dev`). Without it, add
`-DQRAFT_STATIC_EXECUTABLE=OFF -DQRAFT_USE_JEMALLOC=OFF` to the first `cmake`
command. Python 3 is optional: it runs the checks of the counting scheme.

Windows (x86-64): the Windows binary is cross-compiled on Linux with MinGW-w64.
Building on Windows itself is not tested. On Debian 13 or Ubuntu 24.04, in the
cloned repository:

```bash
sudo apt-get update && sudo apt-get install -y g++-mingw-w64-x86-64-posix cmake git
cmake -S . -B build-windows -DCMAKE_SYSTEM_NAME=Windows -DCMAKE_SYSTEM_PROCESSOR=x86_64 \
  -DCMAKE_CXX_COMPILER=x86_64-w64-mingw32-g++-posix -DCMAKE_BUILD_TYPE=Release \
  -DQRAFT_BUILD_TESTS=OFF -DQRAFT_USE_JEMALLOC=OFF
cmake --build build-windows --target qraft -j 4
```

This writes `build-windows/qraft.exe`, a static binary. Copy it to Windows and
run `qraft selftest` there.

macOS is not supported yet.

The release archives are built by `release/build-release.sh` in a Docker image
pinned to a fixed Debian image and package snapshot, so the same commit gives
the same archives:

```bash
release/build-release.sh dist     # needs Docker; JOBS=4 limits the build to 4 CPUs
cd dist && sha256sum -c SHA256SUMS
```

It writes `qraft-1.0.0-linux-x86_64.tar.gz`, `qraft-1.0.0-windows-x86_64.zip`
and `SHA256SUMS`, and checks each binary with its selftest (the Windows one
under Wine). Unpack an archive and run its install check as shown under
[Install](#install).

## Citing QRAFT

Ashrafur Rahman, Anik Saha and Md. Shamsuzzoha Bayzid. QRAFT: Scaling
coalescent-based species tree inference to one million taxa. Manuscript in
preparation.

`CITATION.cff` holds the same reference in machine-readable form. If you
simulate the datasets, please also cite SimPhy-X (Anik Saha, Ashrafur Rahman and
Md. Shamsuzzoha Bayzid, in preparation; <https://github.com/aaniksahaa/simphy-x>)
and SimPhy (Diego Mallo, Leonardo de Oliveira Martins and David Posada, SimPhy:
phylogenomic simulation of gene, locus, and species trees, Systematic Biology
65(2):334-344, 2016).

## Data

Outputs and logs of the paper's experiments:
<https://huggingface.co/datasets/risenfromashes/qraft-data> (CC BY 4.0), in
preparation. The simulated gene trees can be regenerated with the scripts in
`simulation/`; see [Simulating the datasets](#simulating-the-datasets).

## License

Copyright (C) 2026 Ashrafur Rahman, Anik Saha and Md. Shamsuzzoha Bayzid.

QRAFT is free software under the GNU Affero General Public License, version 3
or (at your option) any later version; see `LICENSE`. The release binary
includes third-party code under its own licenses; see
`THIRD_PARTY_NOTICES.md`.
