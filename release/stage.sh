# Steps shared by the release scripts (sourced, after names.sh).

# release_version SRC prints the version in SRC/CMakeLists.txt.
release_version() {
  sed -n 's/^project(qraft VERSION \([0-9.]*\).*/\1/p' "$1/CMakeLists.txt"
}

# stage_common DIR SRC copies the files every archive holds into DIR: the
# license, the documents, the example and Taskflow's license text.
stage_common() {
  local dir=$1 src=$2
  mkdir -p "$dir/examples" "$dir/licenses"
  cp "$src/LICENSE" "$src/README.md" "$src/THIRD_PARTY_NOTICES.md" \
     "$src/CITATION.cff" "$src/release/RELEASE_NOTES.md" "$dir/"
  cp "$src"/examples/* "$dir/examples/"
  cp "$src/external/taskflow/LICENSE" "$dir/licenses/taskflow.txt"
}

# stage_unix DIR QRAFT SRC lays out a Linux or macOS archive in DIR: the
# binary QRAFT as bin/qraft, its man page and shell completions, the common
# files and install-check.sh.
stage_unix() {
  local dir=$1 qraft=$2 src=$3 shell
  mkdir -p "$dir/bin" "$dir/man/man1" "$dir/completions"
  cp "$qraft" "$dir/bin/qraft"
  "$qraft" help man > "$dir/man/man1/qraft.1"
  for shell in bash zsh fish; do
    "$qraft" completion "$shell" > "$dir/completions/qraft.$shell"
  done
  stage_common "$dir" "$src"
  cp "$src/release/install-check.sh" "$dir/"
}

# pack TARGET VERSION STAGE OUT SRC writes the archive of TARGET to OUT from
# its staged folder in STAGE (release/make_archive.py).
pack() {
  local target=$1 version=$2 stage=$3 out=$4 src=$5
  python3 "$src/release/make_archive.py" \
    "$out/$(archive_file "$target" "$version")" "$stage" \
    "$(archive_name "$target" "$version")"
}
