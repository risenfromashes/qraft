#!/usr/bin/env bash
# Writes OUT: the SHA-256 of the trees and branch supports one qraft binary
# infers on the test fixtures, one line per run. The release workflow writes
# it for each of the four release binaries and requires the four files to be
# identical; tests/cli_ablations.cmake pins the trees of every mode on its own.
# Run from the repository root.
#
#   release/tree-hashes.sh QRAFT OUT
set -euo pipefail
qraft=$1
out=$2
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
sha256() {
  if command -v sha256sum > /dev/null; then sha256sum "$1"; else shasum -a 256 "$1"; fi |
    cut -d ' ' -f 1
}
# run NAME ARGS... runs qraft ARGS... -o TREE and prints "HASH  NAME".
run() {
  local name=$1
  shift
  "$qraft" "$@" -o "$tmp/$name.tre" -t 4 -q
  printf '%s  %s\n' "$(sha256 "$tmp/$name.tre")" "$name"
}
complete=tests/trees/astral2_200taxa_50genes.tre
many=tests/trees/astral2_200taxa_1000genes.tre
missing=tests/trees/astral2_200taxa_50genes_25missing.tre
example=examples/genes-100taxa.tre
{
  run example -i "$example"
  run example-support -i "$example" --support localpp
  run genes50 -i "$complete"
  run genes50-extended -i "$complete" --extended
  run genes50-support -i "$complete" --support localpp
  run genes50-support-detail -i "$complete" --support localpp --support-detail \
    --support-frequencies quartet
  run genes1000 -i "$many"
  run genes1000-support -i "$many" --support localpp
  run genes1000-crossover -i "$many" --test-crossover 100
  run missing -i "$missing"
  run missing-crossover -i "$missing" --test-crossover 50
  run support-command support -i "$complete" --tree "$tmp/genes50.tre"
} > "$out"
cat "$out"
