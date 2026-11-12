#!/usr/bin/env bash
# Writes DIR/SHA256SUMS: one line per release archive in DIR (.tar.gz and
# .zip), in name order, in the format `sha256sum -c` reads.
#
#   release/checksums.sh DIR
set -euo pipefail
cd "$1"
shopt -s nullglob
archives=(*.tar.gz *.zip)
if [ ${#archives[@]} -eq 0 ]; then
  echo "checksums.sh: no archives in $1" >&2
  exit 1
fi
printf '%s\n' "${archives[@]}" | LC_ALL=C sort | xargs sha256sum > SHA256SUMS
cat SHA256SUMS
