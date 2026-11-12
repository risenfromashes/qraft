# QRAFT 1.0.1

- `--guide-tree` now always seeds the starting tree: the majority-rule
  extended consensus leaves out the gene-tree clusters that conflict with the
  guide, and the guide's bipartitions then refine its polytomies. In 1.0.0
  this needed `--guide-seed`; without it the guide edited the consensus after
  it was built, which no longer happens. Runs that used `--guide-tree
  --guide-seed` give the same trees as before; runs that used `--guide-tree`
  alone now give the seeded trees.
- `--guide-seed` has no effect any more and prints a note; it is still
  accepted so that older command lines run.
- `--guide-tree` with `--starting-tree` is now a usage error (exit status 2):
  the guide seeds QRAFT's own consensus, so give one or the other.

## Release archives

- `qraft-1.0.1-linux-x86_64.tar.gz`: a static binary for x86-64 Linux, the
  man page, shell completions, a small example and `install-check.sh`.
- `qraft-1.0.1-windows-x86_64.zip`: a static `qraft.exe` for x86-64 Windows,
  the same example and `install-check.cmd`.

Each holds the license, the third-party notices with the license texts of the
code linked into its binary (`licenses/`) and `BUILDINFO.txt`, which names the
compiler and packages. `SHA256SUMS` holds the checksums of both. The two
binaries infer the same trees.

```bash
sha256sum -c SHA256SUMS
tar xf qraft-1.0.1-linux-x86_64.tar.gz
qraft-1.0.1-linux-x86_64/install-check.sh
```

# QRAFT 1.0.0

The first release.

- Species-tree inference from gene trees by quartet regrafting: the
  majority-rule extended consensus, a recursive divide with a statistical
  test of each bipartition, contextual regrafting and NNI refinement, with
  every quartet counted exactly. The same tree for any number of threads.
- Extended mode (`--mode extended`): unbounded batched NNI and an exact SPR
  refinement, for noisy biological data.
- Guide trees (`--guide-tree`), branch support (`--support localpp`, and
  `qraft support` for a given tree), gene trees that miss taxa or repeat them.
- The ablations of the paper as options, among them the refined MRE
  (`--refine-start`), which gives in one run the trees of the four-run chain
  that `qraft help ablations` lists.
- Commands `check`, `version`, `selftest`, `completion` and `help`; input
  errors are reported before the run; `-q` and `-v`; exit status 0, 1, 2 or 3.

## Release archives

- `qraft-1.0.0-linux-x86_64.tar.gz`: a static binary for x86-64 Linux, the
  man page, shell completions, a small example and `install-check.sh`.
- `qraft-1.0.0-windows-x86_64.zip`: a static `qraft.exe` for x86-64 Windows,
  the same example and `install-check.cmd`.

Each holds the license, the third-party notices with the license texts of the
code linked into its binary (`licenses/`) and `BUILDINFO.txt`, which names the
compiler and packages. `SHA256SUMS` holds the checksums of both. The two
binaries infer the same trees.

```bash
sha256sum -c SHA256SUMS
tar xf qraft-1.0.0-linux-x86_64.tar.gz
qraft-1.0.0-linux-x86_64/install-check.sh
```
