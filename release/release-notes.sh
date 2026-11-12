#!/usr/bin/env bash
# Prints the section of release/RELEASE_NOTES.md for VERSION, without its
# title line, for the release page; fails if there is none.
#
#   release/release-notes.sh VERSION
set -euo pipefail
notes=$(awk -v title="# QRAFT $1" '
  $0 == title { found = 1; next }
  found && /^# QRAFT / { exit }
  found { print }
' "$(dirname "$0")/RELEASE_NOTES.md")
if [ -z "${notes//[[:space:]]/}" ]; then
  echo "release-notes.sh: no section '# QRAFT $1' in RELEASE_NOTES.md" >&2
  exit 1
fi
# Without the blank lines at its start.
printf '%s\n' "$notes" | sed -e '/[^[:space:]]/,$!d'
