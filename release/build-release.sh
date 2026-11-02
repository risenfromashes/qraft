#!/usr/bin/env bash
# Builds the static Linux x86-64 release archive in a pinned Docker image.
#
#   release/build-release.sh [OUT_DIR]
#
# The archive and SHA256SUMS go to OUT_DIR (default dist/). Times inside the
# archive come from the commit date (SOURCE_DATE_EPOCH), so two builds of the
# same commit with the same image give the same archive. JOBS=N limits the
# build to N compile jobs and N CPUs (default: every CPU).
set -euo pipefail
root=$(git -C "$(dirname "$0")/.." rev-parse --show-toplevel)
out=${1:-$root/dist}
mkdir -p "$out"
# Docker mounts only absolute paths; a relative one would name a volume.
out=$(cd "$out" && pwd)
if [ -n "$(git -C "$root" status --porcelain --untracked-files=no)" ]; then
  echo "warning: the working tree has uncommitted changes; the version will say -dirty" >&2
fi
git -C "$root" submodule update --init external/taskflow
epoch=$(git -C "$root" log -1 --format=%ct)
docker build -t qraft-release-build "$root/release"
limits=()
if [ -n "${JOBS:-}" ]; then
  limits=(--cpus "$JOBS")
fi
docker run --rm "${limits[@]}" \
  --user "$(id -u):$(id -g)" \
  -e SOURCE_DATE_EPOCH="$epoch" -e HOME=/tmp -e JOBS="${JOBS:-}" \
  -v "$root:/repo:ro" -v "$out:/out" \
  qraft-release-build bash /repo/release/inside.sh
echo "release files in $out:"
ls -l "$out"
