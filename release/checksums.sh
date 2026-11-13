#!/usr/bin/env bash
# Writes DIR/SHA256SUMS: one line per release file in DIR (.tar.gz, .zip, .deb
# and .rpm), in name order, in the format `sha256sum -c` reads.
#
#   release/checksums.sh DIR
set -euo pipefail
cd "$1"
shopt -s nullglob
files=(*.tar.gz *.zip *.deb *.rpm)
if [ ${#files[@]} -eq 0 ]; then
  echo "checksums.sh: no release files in $1" >&2
  exit 1
fi
printf '%s\n' "${files[@]}" | LC_ALL=C sort | xargs sha256sum > SHA256SUMS
cat SHA256SUMS
