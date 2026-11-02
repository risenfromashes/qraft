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
