#!/bin/sh
# Checks an unpacked QRAFT release: the built-in selftest, then the example
# dataset against the tree this release infers from it.
set -eu
here=$(cd "$(dirname "$0")" && pwd)
qraft="$here/bin/qraft"
"$qraft" --version
"$qraft" selftest
out=$(mktemp)
trap 'rm -f "$out"' EXIT
"$qraft" -i "$here/examples/genes-100taxa.tre" -o "$out" -q
if cmp -s "$out" "$here/examples/expected-100taxa.tre"; then
  echo "example: same tree as expected"
  echo "install check passed"
else
  echo "example: the tree differs from examples/expected-100taxa.tre" >&2
  echo "install check FAILED" >&2
  exit 1
fi
